"Somewhere where cells can live. This interface assumes some fields are implemented."
abstract type Environment{T <: AbstractCell} end
"Returns a vector of all cells in the environment (both dead and alive)."
getcells(env::Environment) = env.cells
getcell(env::Environment, sigma) = getcells(env)[sigma]
"Returns an iterable of cells living in the environment."
livingcells(env::Environment) = (cell for cell in getcells(env) if isalive(cell))
getmatrix(env::Environment) = env.matrix
getedgeset(env::Environment) = env.edgeset
getsigma(env::Environment, pos::MatrixPos) = getmatrix(env)[pos]
setsigma(env::Environment, pos::MatrixPos, val) = getmatrix(env)[pos] = val
nextsigma(env::Environment) = lastindex(getcells(env)) + 1
celltype(::Environment{T})  where T = T

"A simple petri dish with not much going on."
struct Dish{T} <: Environment{T}
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

function spawncell!(env::Environment,
                    tau::Int, 
                    cell_len::Integer, 
                    divtime::Integer, 
                    spawncenter=getrandompos(getmatrix(env), ceil(Int, cell_len / 2)))
    sigma = nextsigma(env)
    bb, area, center = initcellpositions!(getmatrix(env), sigma, spawncenter, cell_len)
    push!(getcells(env), celltype(env)(sigma, tau, area, center, divtime))
    union!(getedgeset(env), getedges(env, bb))
end

"Initialize a cell on a matrix and return its area and bounding box."
function initcellpositions!(matrix::Matrix, sigma::Integer, center::MatrixPos, cell_len::Integer)
    bb = BoundingBox(center, center)
    # Recalculate area and center to account for positions that couldnt be initialized
    area = 0
    centerx, centery = 0., 0.
    rangebb = BoundingBox(
        MatrixPos(getx(center) - floor(cell_len / 2), 
                  gety(center) - floor(cell_len / 2)),
        MatrixPos(getx(center) + ceil(cell_len / 2) - 1, 
                  gety(center) + ceil(cell_len / 2) - 1)
    )
    for pos in iterpositions(rangebb)
        if matrix[pos] == 0
            matrix[pos] = sigma
            area += 1
            centerx += getx(pos) / area
            centery += gety(pos) / area
            addpos!(bb, pos)
        end
    end
    bb, area, Pos(centerx, centery)
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

function updateedges!(env::Environment, newsigma, pos)
    for neigh in moore_neighbors(pos)
        sigmaneigh = getsigma(env, neigh)
        if newsigma == sigmaneigh
            removeedges!(getedgeset(env), pos, neigh)
        elseif sigmaneigh >= 0
            addedges!(getedgeset(env), pos, neigh)
            addedges!(getedgeset(env), pos, neigh)
    end end
end

function updateattributes!(env::Environment, oldsigma, newsigma, edge)
    if newsigma > 0
        newcell = getcell(env, newsigma)
        addarea!(newcell, 1)
        addmomentum!(newcell, edge[2])
    end
    if oldsigma > 0
        oldcell = getcell(env, oldsigma)
        addarea!(oldcell, -1)
        removemomentum!(oldcell, edge[2])
    end
end