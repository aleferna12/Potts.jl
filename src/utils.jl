"A position on the simulation matrix."
struct Pos{T<:Number}
    x::T
    y::T
end
inbounds(pos::Pos, bounds::Tuple) = 0 < pos.x <= bounds[1] && 0 < pos.y <= bounds[2]
filterinbounds(positions, bounds::Tuple) = filter(pos -> inbounds(pos, bounds), positions)

"Often in CPM simulations positions are represented as matrix indices, so we provide this alias definition."
MatrixPos = Pos{Int}
Base.getindex(array::Array, pos::MatrixPos) = array[pos.x, pos.y]
Base.setindex!(array::Array, value, pos::MatrixPos) = array[pos.x, pos.y] = value
Base.checkindex(::Type{Bool}, inds::AbstractUnitRange, index::MatrixPos) = 0 < index.x * index.y <= length(inds)
getadjacentx(pos::MatrixPos) = [MatrixPos(pos.x - 1, pos.y), MatrixPos(pos.x + 1, pos.y)]
getadjacenty(pos::MatrixPos) = [MatrixPos(pos.x, pos.y - 1), MatrixPos(pos.x, pos.y + 1)]
vonneumann_neighbors(pos::MatrixPos) = [getadjacentx(pos); getadjacenty(pos)]
moore_neighbors(pos::MatrixPos) = [vonneumann_neighbors(pos)...,
                                   MatrixPos(pos.x - 1, pos.y - 1),
                                   MatrixPos(pos.x - 1, pos.y + 1),
                                   MatrixPos(pos.x + 1, pos.y - 1),
                                   MatrixPos(pos.x + 1, pos.y + 1)]

function getrandompos(matrix::Matrix, borderpadding::Integer=0)
    xrange = range(1 + borderpadding, size(matrix)[1] - borderpadding)
    yrange = range(1 + borderpadding, size(matrix)[2] - borderpadding)
    MatrixPos(rand(xrange), rand(yrange))
end

Edge = Pair{MatrixPos, MatrixPos}

function removeedges!(edgeset, pos1, pos2) 
    delete!(edgeset, Edge(pos1, pos2))
    delete!(edgeset, Edge(pos2, pos1))
end

function addedges!(edgeset, pos1, pos2)
    push!(edgeset, Edge(pos1, pos2))
    push!(edgeset, Edge(pos2, pos1))
end

allsame(x) = all(y -> y == first(x), x)

orderpair(x1, x2) = x2 < x1 ? x2 => x1 : x1 => x2