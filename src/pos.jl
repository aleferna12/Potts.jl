"A position on the simulation matrix."
struct Pos{T<:Number}
    x::T
    y::T
end
getx(pos::Pos) = pos.x
gety(pos::Pos) = pos.y
Base.round(postype::Type, pos::Pos) = postype(round(getx(pos)), round(gety(pos)))
Base.round(pos::Pos) = round(typeof(pos), pos)
Base.convert(postype::Type{<:Pos}, pos::NTuple{2}) = postype(pos[1], pos[2])

"Determines on which side of a line (specified by 'm' and 'n') a point lies. Returns 1 for 'up' or 'left' and -1 otherwise."
function whichside(m, n, pos::Pos)
    y = m * getx(pos) + n
    sign(gety(pos) - y)
end

"Often in CPM simulations positions are represented as matrix indices, so we provide this alias definition."
const MatrixPos = Pos{Int}
Base.getindex(matrix::Matrix, pos::MatrixPos) = matrix[pos.x, pos.y]
Base.getindex(matrix::Matrix, positions::MatrixPos...) = [matrix[pos] for pos in positions]
Base.setindex!(matrix::Matrix, value, pos::MatrixPos) = matrix[pos.x, pos.y] = value
Base.setindex!(matrix::Matrix, value, positions::MatrixPos...) = for pos in positions matrix[pos] = value end
Base.checkbounds(::Type{Bool}, matrix::Matrix, index::MatrixPos) = checkbounds(Bool, matrix, CartesianIndex(getx(index), gety(index)))
Base.convert(::Type{MatrixPos}, pos::Pos) = MatrixPos(getx(pos), gety(pos))
moore_neighbors(pos::MatrixPos) = [MatrixPos(pos.x - 1, pos.y),
                                   MatrixPos(pos.x + 1, pos.y),
                                   MatrixPos(pos.x, pos.y - 1),
                                   MatrixPos(pos.x, pos.y + 1),
                                   MatrixPos(pos.x - 1, pos.y - 1),
                                   MatrixPos(pos.x - 1, pos.y + 1),
                                   MatrixPos(pos.x + 1, pos.y - 1),
                                   MatrixPos(pos.x + 1, pos.y + 1)]

function getrandompos(matrix::Matrix, borderpadding::Integer=0)
    xrange = range(1 + borderpadding, size(matrix)[1] - borderpadding)
    yrange = range(1 + borderpadding, size(matrix)[2] - borderpadding)
    MatrixPos(rand(xrange), rand(yrange))
end

const Edge = Pair{MatrixPos, MatrixPos}

function removeedges!(edgeset, pos1, pos2)
    pop!(edgeset, Edge(pos1, pos2), Edge(pos2, pos1))
end

function addedges!(edgeset, pos1, pos2)
    push!(edgeset, Edge(pos1, pos2), Edge(pos2, pos1))
end