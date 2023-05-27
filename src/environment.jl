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

function updateattributes!(env::Environment, oldsigma, newsigma, edge)
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