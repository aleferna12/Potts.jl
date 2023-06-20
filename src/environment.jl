"Somewhere where cells can live."
abstract type Environment{T <: AbstractCell} end

function nextsigma(env::Environment)
    for (i, cell) in enumerate(getcells(env))  
        if !isalive(cell)
            return i
    end end  
    lastindex(getcells(env)) + 1
end

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

function addcell!(env::Environment, cell::AbstractCell)
    if getsigma(cell) > length(getcells(env))
        push!(getcells(env), cell)
    else
        getcells(env)[getsigma(cell)] = cell
    end
    nothing
end

"Spawn a cell on a set of specified positions."
function spawncell!(env::Environment{CT}, positions; cellkwargs...) where {CT <: AbstractCell}
    sigma = nextsigma(env)
    area = length(positions)
    centerx = mean(getx(pos) for pos in positions)
    centery = mean(gety(pos) for pos in positions)
    newcell = CT(sigma=sigma, area=area, center=(centerx, centery); cellkwargs...)
    addcell!(env, newcell)
    for pos in positions
        setsigma!(env, pos, sigma)
    end
    union!(getedgeset(env), getedges(env, containing(positions)))
    newcell
end

"Initializes a square-shaped cell on free positions on the matrix. This is used during the model setup."
function initcell!(env::Environment,
                   cell_len::Integer,
                   spawncenter::MatrixPos;
                   cellkwargs...)
    bb = BoundingBox(MatrixPos(getx(spawncenter) - floor(cell_len / 2), 
                               gety(spawncenter) - floor(cell_len / 2)),
                     MatrixPos(getx(spawncenter) + ceil(cell_len / 2) - 1, 
                               gety(spawncenter) + ceil(cell_len / 2) - 1))
    inbounds = (pos for pos in iterpositions(bb) if checkbounds(Bool, getmatrix(env), pos))
    positions = collect(emptypositions(getmatrix(env), inbounds))
    if !isempty(positions)
        spawncell!(env, positions; cellkwargs...)
end end

function emptypositions(matrix::Matrix, positions)
    (pos for pos in positions if matrix[pos] == MED)
end

function getedges(env::Environment, bb::BoundingBox)
    edges = Set{Edge}()
    for pos in iterpositions(bb)
        sigmapos = getsigma(env, pos)
        for neigh in moore_neighbors(pos)
            sigmaneigh = getsigma(env, neigh)
            if sigmapos != BORDER && sigmaneigh != BORDER && sigmapos != sigmaneigh
                addedges!(edges, pos, neigh)
    end end end
    edges
end

function iterpositions(env::Environment, cell::AbstractCell)
    (pos for pos in CellPosIter(getsigma(cell), getarea(cell), round(getcenter(cell)), getmatrix(env)))
end

function updateattributes!(env::Environment, targetcellarea, divtargetcellarea)
    for cell in livingcells(env)
        ta, divta = targetcellarea[gettau(cell)], divtargetcellarea[gettau(cell)]
        settargetarea!(cell, growthperc(cell) * (divta - ta) + ta)
end end

function updateattributes!(env::Environment{<:EvolvableCell}, targetcellarea, divtargetcellarea)
    @invoke updateattributes!(env::Environment, targetcellarea, divtargetcellarea)
    for cell in livingcells(env)
        update!(getgenome(cell))
end end

function select!(env::Environment{CT}) where CT
    killed = CT[]
    for cell in livingcells(env)
        if getarea(cell) <= 0
            kill!(env, cell)
            push!(killed, cell)
    end end
    killed
end

function kill!(env::Environment, cell::AbstractCell)
    if getarea(cell) > 0
        for pos in iterpositions(env, cell)
            setsigma!(env, pos, MED)
    end end
    kill!(cell)
end

"Checks which cells can reproduce and divide them."
function reproduce!(env::Environment{CT}; divkwargs...) where CT
    mothers, daughters = CT[], CT[]
    for cell in livingcells(env)
        if isdividing(cell) && dividenow!(cell)
            cells = dividecell!(env, cell; divkwargs...)
            if !isnothing(cells)
                push!(mothers, cells[1])
                push!(daughters, cells[2])
    end end end
    mothers, daughters
end

"Tells a cell to divide in two. If the cell is too small this function won't do anything."
function dividecell!(env::Environment, cell::AbstractCell)
    oldpositions, newpositions = split(iterpositions(env, cell))

    if !isempty(oldpositions) && !isempty(newpositions)
        newcell = spawncell!(env, newpositions, tau=gettau(cell), targetarea=gettargetarea(cell), divtimer=IterationTimer(getdivtime(cell)))

        # Update old cell's attributes
        for pos in newpositions
            losepos!(cell, pos)
        end
        return cell, newcell
end end

function dividecell!(env::Environment, cell::EvolvableCell; mu, mustd)
    oldpositions, newpositions = split(iterpositions(env, cell))

    if !isempty(oldpositions) && !isempty(newpositions)
        newcell = spawncell!(env, newpositions, tau=gettau(cell), targetarea=gettargetarea(cell), divtimer=IterationTimer(getdivtime(cell)), genome=deepcopy(getgenome(cell)))
        mutate!(getgenome(newcell), mu, mustd)

        # Update old cell's attributes
        for pos in newpositions
            losepos!(cell, pos)
        end
        return cell, newcell
end end

"A simple petri dish with not much going on."
struct Dish{T} <: Environment{T}
    cells::Vector{T}
    matrix::Matrix{Int}
    edgeset::Set{Edge}
end
Dish(celltype::Type{T}, fieldsize::Integer) where T <: AbstractCell = Dish(celltype[], createcellmatrix(fieldsize), Set{Edge}())
"Returns a vector of all cells in the environment (both dead and alive)."
getcells(dish::Dish) = dish.cells
getcell(dish::Dish, sigma) = getcells(dish)[sigma]
getmatrix(dish::Dish) = dish.matrix
"Returns an iterable of cells living in the environment."
livingcells(dish::Dish) = (cell for cell in getcells(dish) if isalive(cell))
getedgeset(dish::Dish) = dish.edgeset
getsigma(dish::Dish, pos::MatrixPos) = getmatrix(dish)[pos]
setsigma!(dish::Dish, pos::MatrixPos, val) = getmatrix(dish)[pos] = val

struct CellPosIter
    sigma::Int
    area::Int
    center::MatrixPos
    matrix::Matrix{Int}
end
Base.length(cpi::CellPosIter) = cpi.area

# Attempt at an efficient way of querying bigger and bigger rings of positions until all cell positions are found
# Maybe an improvement would be to first iterate over all positions in a bb that likely contains all cell positions and only then expand on that
function Base.iterate(cpi::CellPosIter, state=(found=0, curopi=OuterPosIter(cpi.center, 2, 2), opistate=(curpos=cpi.center, dir=:up)))
    found, curopi, opistate = state.found, state.curopi, state.opistate
    if found == cpi.area
        return nothing
    end

    while (it = iterate(curopi, opistate)) !== nothing
        if checkbounds(Bool, cpi.matrix, it[1]) && cpi.matrix[it[1]] == cpi.sigma
            return it[1], (found=found + 1, curopi=curopi, opistate=it[2])
        end
        opistate = it[2]
    end
    minpos = MatrixPos(getx(curopi.minpos) - 1, gety(curopi.minpos) - 1)
    iterate(cpi, (found=found, curopi=OuterPosIter(minpos, curopi.lengthx + 2, curopi.lengthy + 2), opistate=(curpos=minpos, dir=:up)))
end