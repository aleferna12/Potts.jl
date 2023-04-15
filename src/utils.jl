struct Pos{T<:Number}
    x::T
    y::T
end
Base.:(<)(pos1::Pos, pos2::Pos) = pos1.x < pos2.x || (pos1.x == pos2.x && pos1.y < pos2.y)

MatrixPos = Pos{Int}
getadjacentx(pos::MatrixPos) = MatrixPos(pos.x - 1, pos.y), MatrixPos(pos.x + 1, pos.y)
getadjacenty(pos::MatrixPos) = MatrixPos(pos.x, pos.y - 1), MatrixPos(pos.x, pos.y + 1)
vonneumann_neighbors(pos::MatrixPos) = getadjacentx(pos)..., getadjacenty(pos)...
moore_neighbors(pos::MatrixPos) = vonneumann_neighbors(pos)...,
                                  MatrixPos(pos.x - 1, pos.y - 1),
                                  MatrixPos(pos.x - 1, pos.y + 1),
                                  MatrixPos(pos.x + 1, pos.y - 1),
                                  MatrixPos(pos.x + 1, pos.y + 1)

function getrandompos(matrix::Matrix, borderpadding::Integer=0)
    xrange = range(1 + borderpadding, size(matrix)[1] - borderpadding)
    yrange = range(1 + borderpadding, size(matrix)[2] - borderpadding)
    MatrixPos(rand(xrange), rand(yrange))
end

struct Edge
    pos1::MatrixPos
    pos2::MatrixPos
    Edge(pos1::MatrixPos, pos2::MatrixPos) = pos2 < pos1 ? new(pos2, pos1) : new(pos1, pos2)
end