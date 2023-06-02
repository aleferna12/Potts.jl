mutable struct BoundingBox
    minpos::MatrixPos
    maxpos::MatrixPos
    BoundingBox(minpos::MatrixPos, maxpos::MatrixPos) = getx(maxpos) >= getx(minpos) && gety(maxpos) >= gety(minpos) ? 
                                                        new(minpos, maxpos) :
                                                        error("negative area in BoundingBox")
end
BoundingBox(minpos, maxpos) = BoundingBox(convert(MatrixPos, minpos), convert(MatrixPos, maxpos))
"Creates a bounding box around a single position."
BoundingBox(pos::MatrixPos) = BoundingBox(pos, pos)
getminpos(bb::BoundingBox) = bb.minpos
getmaxpos(bb::BoundingBox) = bb.maxpos
getsidex(bb::BoundingBox) = getx(bb.maxpos) - getx(bb.minpos) + 1
getsidey(bb::BoundingBox) = gety(bb.maxpos) - gety(bb.minpos) + 1
getarea(bb::BoundingBox) = getsidex(bb) * getsidey(bb)
getvertices(bb::BoundingBox) = [bb.minpos, MatrixPos(bb.minpos.x, bb.maxpos.y), bb.maxpos, MatrixPos(bb.maxpos.x, bb.minpos.y)]
getcenter(bb::BoundingBox) = Pos((bb.minpos.x + bb.maxpos.x) / 2, (bb.minpos.y + bb.maxpos.y) / 2)
"Iterate over every position inside of the BoundingBox."
iterpositions(bb::BoundingBox) = (MatrixPos(x, y) for x in bb.minpos.x:bb.maxpos.x for y in bb.minpos.y:bb.maxpos.y)
function iterouterpositions(bb::BoundingBox)
    if getsidex(bb) < 3 || getsidey(bb) < 3
        return iterpositions(bb)
    end
    (pos for pos in OuterPosIter(getminpos(bb), getsidex(bb), getsidey(bb)))
end

"Creates a BoundingBox that contains a set of positions."
function containing(positions)
    bb = BoundingBox(first(positions))
    for pos in positions
        addpos!(bb, pos)
    end
    bb
end

"Add a position to a BoundingBox by enlarging it."
function addpos!(bb::BoundingBox, pos::MatrixPos)
    bb.minpos = Pos(min(bb.minpos.x, pos.x), min(bb.minpos.y, pos.y))
    bb.maxpos = Pos(max(bb.maxpos.x, pos.x), max(bb.maxpos.y, pos.y))
end

struct OuterPosIter
    minpos::MatrixPos
    lengthx::Int
    lengthy::Int
    OuterPosIter(minpos, lengthx, lengthy) = lengthx > 1 && lengthy > 1 ? new(minpos, lengthx, lengthy) : error("are smaller than 4")
end
Base.length(opi::OuterPosIter) = 2 * (opi.lengthx + opi.lengthy) - 4

# Somewhat an ugly solution but is orders of magnitude faster than using Iterators.flatten on four iters
function Base.iterate(opi::OuterPosIter, state=(curpos=opi.minpos, dir=:up))
    dir, curpos = state.dir, state.curpos
    if dir === :left
        if getx(curpos) - getx(opi.minpos) <= 0
            return nothing
        end
        return curpos, (curpos=MatrixPos(getx(curpos) - 1, gety(curpos)), dir=:left)
    end
    if dir === :up
        if gety(curpos) - gety(opi.minpos) >= opi.lengthy - 1
            return curpos, (curpos=MatrixPos(getx(curpos) + 1, gety(curpos)), dir=:right)
        end
        return curpos, (curpos=MatrixPos(getx(curpos), gety(curpos) + 1), dir=:up)
    end
    if dir === :right
        if getx(curpos) - getx(opi.minpos) >= opi.lengthx - 1
            return curpos, (curpos=MatrixPos(getx(curpos), gety(curpos) - 1), dir=:down)
        end
        return curpos, (curpos=MatrixPos(getx(curpos) + 1, gety(curpos)), dir=:right)
    end
    if dir === :down
        if gety(curpos) - gety(opi.minpos) <= 0
            return curpos, (curpos=MatrixPos(getx(curpos) - 1, gety(curpos)), dir=:left)
        end
        return curpos, (curpos=MatrixPos(getx(curpos), gety(curpos) - 1), dir=:down)
end end