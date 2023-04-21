@enum Tau veg=1 div

mutable struct AttrSet
    sigma::Int
    tau::Tau
    area::Int
    center::Pos
end

mutable struct Genome{}
    sigma::Int
end

struct Cells
    attrsets::Vector{AttrSet}
    genomes::Vector{Genome}
    matrix::Matrix{Int}
    edgeset::Set{Edge}
end
Cells(size::Integer) = Cells(AttrSet[], Genome[], zeros(Int, size, size), Set{Edge}())
Base.length(cells::Cells) = length(cells.attrsets)
Base.lastindex(cells::Cells) = lastindex(cells.attrsets)

function spawncell!(cells::Cells, cell_len::Integer)
    sigma = lastindex(cells) + 1
    spawncenter = getrandompos(cells.matrix, cell_len)
    area, bb = initcellpositions!(cells.matrix, sigma, spawncenter, cell_len)
    attrset = AttrSet(
        sigma,
        veg,
        area,
        getcenter(bb)
    )
    push!(cells.attrsets, attrset)
    union!(cells.edgeset, getedges(bb))
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