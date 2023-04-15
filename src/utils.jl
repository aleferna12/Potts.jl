struct Pos{T<:Number}
    x::T
    y::T
end
Base.:(<)(pos1::Pos, pos2::Pos) = pos1.x < pos2.x || (pos1.x == pos2.x && pos1.y < pos2.y)
MatrixPos = Pos{Int}

function getrandompos(matrix::Matrix, borderpadding::Integer=0)
    xrange = range(1 + borderpadding, size(matrix)[1] - borderpadding)
    yrange = range(1 + borderpadding, size(matrix)[2] - borderpadding)
    MatrixPos(rand(xrange), rand(yrange))
end

MatrixSite{T<:Number} = Pair{MatrixPos, T}
getadjacent(pos::MatrixPos) = getadjacentx(pos)..., getadjacenty(pos)...
getadjacentx(pos::MatrixPos) = MatrixPos(pos.x - 1, pos.y), MatrixPos(pos.x + 1, pos.y)
getadjacenty(pos::MatrixPos) = MatrixPos(pos.x, pos.y - 1), MatrixPos(pos.x, pos.y + 1)
getmatrixsite(pos::MatrixPos, matrix::Matrix) = pos => matrix[pos.x, pos.y]
getmatrixsite(x::Int, y::Int, matrix::Matrix) = getmatrixsite(MatrixPos(x, y), matrix)

struct Edge{T<:Number}
    site1::MatrixSite{T}
    site2::MatrixSite{T}
    Edge(site1::MatrixSite{S}, site2::MatrixSite{S}) where S<: Number = 
        site2.first < site1.first ? new{S}(site2, site1) : new{S}(site1, site2)
end