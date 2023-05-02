mutable struct BoundingBox
    minpos::MatrixPos
    maxpos::MatrixPos
end
Base.length(bb::BoundingBox) = length(iterpositions(bb))
getvertices(bb::BoundingBox) = [bb.minpos, MatrixPos(bb.minpos.x, bb.maxpos.y), bb.maxpos, MatrixPos(bb.maxpos.x, bb.minpos.y)]
getcenter(bb::BoundingBox) = Pos((bb.minpos.x + bb.maxpos.x) / 2, (bb.minpos.y + bb.maxpos.y) / 2)
"Iterate over every position inside of the BoundingBox."
iterpositions(bb::BoundingBox) = (MatrixPos(x, y) for x in bb.minpos.x:bb.maxpos.x, 
                                                      y in bb.minpos.y:bb.maxpos.y)

"Add a position to a BoundingBox by enlarging it."
function addpos!(bb::BoundingBox, pos::MatrixPos)
    bb.minpos = Pos(min(bb.minpos.x, pos.x), min(bb.minpos.y, pos.y))
    bb.maxpos = Pos(max(bb.maxpos.x, pos.x), max(bb.maxpos.y, pos.y))
end