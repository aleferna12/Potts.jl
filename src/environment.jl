"Somewhere where cells can live."
abstract type Environment end
getcells(env::Environment) = env.cells
getcell(env::Environment, sigma) = getcells(env)[sigma]
getmatrix(env::Environment) = env.matrix
getedgeset(env::Environment) = env.edgeset
getsigma(env::Environment, pos::MatrixPos) = getmatrix(env)[pos]
setsigma(env::Environment, pos::MatrixPos, val) = getmatrix(env)[pos] = val

"A simple petri dish with not much going on."
struct Dish{T <: AbstractCell} <: Environment
    cells::Vector{T}
    matrix::Matrix{Int}
    edgeset::Set{Edge}
end
Dish(celltype::Type{T}, fieldsize::Integer) where T <: AbstractCell = Dish(celltype[], createcellmatrix(fieldsize), Set{Edge}())

function createcellmatrix(fieldsize)
    m = zeros(Int, fieldsize, fieldsize)
    for i in 1:fieldsize  # Borders
        m[i, 1] = -1
        m[i, fieldsize] = -1
        m[1, i] = -1
        m[fieldsize, i] = -1
    end
    m
end

function spawncell!(env::Environment, tau::Int, cell_len::Integer, spawncenter=getrandompos(getmatrix(env), ceil(Int, cell_len / 2)))
    sigma = lastindex(getcells(env)) + 1
    area, bb = initcellpositions!(getmatrix(env), sigma, spawncenter, cell_len)
    push!(getcells(env), Cell(sigma, tau, area, getcenter(bb)))
    union!(getedgeset(env), getedges(env, bb))
end

"Initialize a cell on a matrix and return its area and bounding box."
function initcellpositions!(matrix::Matrix, sigma::Integer, center::MatrixPos, cell_len::Integer)
    area = 0
    bb = BoundingBox(center, center)
    rangebb = BoundingBox(
        MatrixPos(center.x - floor(cell_len / 2), 
                  center.y - floor(cell_len / 2)),
        MatrixPos(center.x + ceil(cell_len / 2) - 1, 
                  center.y + ceil(cell_len / 2) - 1)
    )
    for pos in iterpositions(rangebb)
        if matrix[pos.x, pos.y] == 0
            matrix[pos.x, pos.y] = sigma
            area += 1
            addpos!(bb, pos)
        end
    end
    area, bb
end

function getedges(env::Environment, bb::BoundingBox)
    edges = Set{Edge}()
    for pos in iterpositions(bb)
        sigmapos = getsigma(env, pos)
        for neigh in moore_neighbors(pos)
            sigmaneigh = getsigma(env, neigh)
            if sigmapos > -1 && sigmaneigh > -1 && sigmapos != sigmaneigh
                addedges!(edges, pos, neigh)
    end end end
    edges
end

function setupenv!(env::Environment, params)
    for tau in 1:params.taus
        for _ in 1:params.ncells[tau]
            spawncell!(env, tau, params.cell_length[tau])
    end end
    env
end

function step!(env::Environment, time::Integer, params)
    updatehamiltonian!(env, params)
    # Iterate the matrix and use the sigmas to calculate and update cell attributes
    # updateattributes!(cells)
end

function updatehamiltonian!(env::Environment, params)
    # TODO: Check c++ code to understand how we used to do it
    # I think they are adding to a 'loop' variable but this variable is an int and they are adding floats (so doing nothing essentially)
    tovisit = length(getedgeset(env)) / 8  # 8 because moore_neighbors gives 8 neighbours
    for _ in 1:ceil(tovisit)
        edge = rand(getedgeset(env)) # TODO: this takes very long! Can it be optimized?
        if getsigma(env, edge[1]) == getsigma(env, edge[2])
            continue
        end

        dH = deltahamiltonian(env, edge, params)
        if acceptcopy(dH, params.boltzmanntemp)
            copyspin!(env, edge)
end end end

"Computes the local energy difference in case a lattice copy event (represented by an edge) occurs."
function deltahamiltonian(env::Environment, copyattempt::Edge, params)
    deltaH = 0.
    sigma1, sigma2 = getsigma(env, copyattempt[1]), getsigma(env, copyattempt[2])
    if sigma1 != 0
        deltaH += deltaHsize(getarea(getcell(env, sigma1)), 1, params.targetcellarea[gettau(getcell(env, sigma1))], params.sizelambda)
    end
    if sigma2 != 0
        deltaH += deltaHsize(getarea(getcell(env, sigma2)), -1, params.targetcellarea[gettau(getcell(env, sigma2))], params.sizelambda)
    end
    neighsigmas = [getsigma(env, nsig) for nsig in moore_neighbors(copyattempt[2])]
    deltaH + deltaHadhenergy(env, sigma2, sigma1, neighsigmas, params)
end

deltaHsize(area, deltaarea, targetarea, sizelambda) = sizelambda * deltaarea * (2 * (area - targetarea) + deltaarea)

function deltaHadhenergy(env,
                         oldsigma,
                         newsigma,
                         neighsigmas,
                         params)
    energy = 0.
    for neighsigma in neighsigmas
        energy += getadhenergy(env, newsigma, neighsigma, params)
        energy -= getadhenergy(env, oldsigma, neighsigma, params)
    end
    energy
end

function getadhenergy(env::Environment, sigma1::Integer, sigma2::Integer, params)::Float64
    if sigma1 == sigma2
        return 0.
    end
    sigma1, sigma2 = ordered(sigma1, sigma2)
    if sigma1 < 0
        return sigma2 == 0 ? 0. : params.borderenergy
    elseif sigma1 == 0
        return params.adhesionmedium[gettau(getcell(env, sigma2))]
    end
    params.adhesiontable[gettau(getcell(env, sigma1)), gettau(getcell(env, sigma2))] * 2
end

function acceptcopy(deltaH, bolzmanntemp)
    if deltaH < 0
        return true
    end
    rand() < ℯ^(-deltaH / bolzmanntemp)
end

function copyspin!(env::Environment, edge::Edge)
    sigma1, sigma2 = getsigma(env, edge[1]), getsigma(env, edge[2])
    setsigma(env, edge[2], sigma1)
    updateattributes!(env, sigma2, sigma1, edge)
    updateedges!(env, sigma1, edge[2])
end

function updateedges!(env::Environment, newsigma, pos)
    for neigh in moore_neighbors(pos)
        sigmaneigh = getsigma(env, neigh)
        if newsigma == sigmaneigh
            removeedges!(getedgeset(env), pos, neigh)
        elseif sigmaneigh > -1
            addedges!(getedgeset(env), pos, neigh)
            addedges!(getedgeset(env), pos, neigh)
    end end
end

function updateattributes!(env::Environment, oldsigma, newsigma, edge::Edge)
    if newsigma != 0
        newcell = getcell(env, newsigma)
        addarea!(newcell, 1)
        addmomentum!(newcell, edge[2])
    end
    if oldsigma != 0
        oldcell = getcell(env, oldsigma)
        addarea!(oldcell, -1)
        removemomentum!(oldcell, edge[2])
    end
end