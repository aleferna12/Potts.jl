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

function getedges(bb::BoundingBox, matrix::Matrix)
    edges = Edge[]
    for x in bb.minpos.x:bb.maxpos.x
        site1 = getmatrixsite(x, bb.minpos.y, matrix)
        site2 = getmatrixsite(x, bb.minpos.y - 1, matrix)
        site3 = getmatrixsite(x, bb.maxpos.y, matrix)
        site4 = getmatrixsite(x, bb.maxpos.y + 1, matrix)
        push!(edges, Edge(site1, site2), Edge(site3, site4))
    end
    for y in bb.minpos.y:bb.maxpos.y
        site1 = getmatrixsite(bb.minpos.x, y, matrix)
        site2 = getmatrixsite(bb.minpos.x - 1, y, matrix)
        site3 = getmatrixsite(bb.maxpos.x, y, matrix)
        site4 = getmatrixsite(bb.maxpos.x + 1, y, matrix)
        push!(edges, Edge(site1, site2), Edge(site3, site4))
    end
    push!(
        edges,
        Edge(getmatrixsite(bb.minpos, matrix), getmatrixsite(bb.minpos.x - 1, bb.minpos.y - 1, matrix)),
        Edge(getmatrixsite(bb.maxpos, matrix), getmatrixsite(bb.maxpos.x + 1, bb.maxpos.y + 1, matrix)),
        Edge(getmatrixsite(bb.minpos.x, bb.maxpos.y, matrix), getmatrixsite(bb.minpos.x - 1, bb.maxpos.y + 1, matrix)),
        Edge(getmatrixsite(bb.maxpos.x, bb.minpos.y, matrix), getmatrixsite(bb.maxpos.x + 1, bb.minpos.y - 1, matrix))
    )
    edges
end

"Add a position to a BoundingBox by enlarging it."
function addpos!(bb::BoundingBox, pos::MatrixPos)
    addposx!(bb, pos.x)
    addposy!(bb, pos.y)
end

# When updating the hamiltonian we can speed up the update process by skipping unecessary min and max calls
function addposx!(bb::BoundingBox, x::Integer)
    bb.minpos = Pos(min(bb.minpos.x, x), bb.minpos.y)
    bb.maxpos = Pos(max(bb.maxpos.x, x), bb.maxpos.y)
end

function addposy!(bb::BoundingBox, y::Integer)
    bb.minpos = Pos(bb.minpos.x, min(bb.minpos.y, y))
    bb.maxpos = Pos(bb.maxpos.x, max(bb.maxpos.y, y))
end

"Remove a pos from a BoundingBox."
function removepos!(bb::BoundingBox, pos::MatrixPos)
    removeposx!(bb, pos.x)
    removeposy!(bb, pos.y)
end

function removeposx!(bb::BoundingBox, x::Integer)
end

function removeposy!(bb::BoundingBox, y::Integer)
end