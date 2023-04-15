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

function getedges(bb::BoundingBox)
    edges = Edge[]
    for x in bb.minpos.x:bb.maxpos.x
        push!(
            edges,
            Edge(MatrixPos(x, bb.minpos.y), MatrixPos(x, bb.minpos.y - 1)),
            Edge(MatrixPos(x, bb.maxpos.y), MatrixPos(x, bb.maxpos.y + 1))
        )
    end
    for y in bb.minpos.y:bb.maxpos.y
        push!(
            edges,
            Edge(MatrixPos(bb.minpos.x, y), MatrixPos(bb.minpos.x - 1, y)),
            Edge(MatrixPos(bb.maxpos.x, y), MatrixPos(bb.maxpos.x + 1, y))
        )
    end
    push!(
        edges,
        Edge(bb.minpos, MatrixPos(bb.minpos.x - 1, bb.minpos.y - 1)),
        Edge(bb.maxpos, MatrixPos(bb.maxpos.x + 1, bb.maxpos.y + 1)),
        Edge(MatrixPos(bb.minpos.x, bb.maxpos.y), MatrixPos(bb.minpos.x - 1, bb.maxpos.y + 1)),
        Edge(MatrixPos(bb.maxpos.x, bb.minpos.y), MatrixPos(bb.maxpos.x + 1, bb.minpos.y - 1))
    )
    edges
end

# Note that to remove a position we would need to check the values from _posx and _posy in the associated BoundaryMap
"Add a position to a BoundingBox by enlarging it."
function addpos!(bb::BoundingBox, pos::MatrixPos)
    bb.minpos = Pos(min(bb.minpos.x, pos.x), min(bb.minpos.y, pos.y))
    bb.maxpos = Pos(max(bb.maxpos.x, pos.x), max(bb.maxpos.y, pos.y))
end

# TODO do I need this? Ideally I wouldn't need to iterate over specific cell lattice sites for anything
# In the C++ code this kind of iteration was used to compute the vectors of motion and the food consumption
struct BoundaryMap
    bb::BoundingBox
    _posx::Vector{Int} # Keeps all the interface positions with the medium and other cells
    _posy::Vector{Int}
end