using Statistics: mean

function initcellpositions!(cellmatrix::Matrix, sigma::Integer, center::MatrixPos, cell_len::Integer)
    area = 0
    bb = BoundingBox(center, center)
    rangebb = BoundingBox(
        MatrixPos(center.x - floor(cell_len / 2), 
                  center.y - floor(cell_len / 2)),
        MatrixPos(center.x + ceil(cell_len / 2) - 1, 
                  center.y + ceil(cell_len / 2) - 1)
    )
    for pos in iterpositions(rangebb)
        if cellmatrix[pos.x, pos.y] == 0
            cellmatrix[pos.x, pos.y] = sigma
            area += 1
            addpos!(bb, pos)
        end
    end
    area, bb
end

function spawncell!(cells::Cells, cellmatrix::Matrix, cell_len::Integer)
    sigma = lastindex(cells) + 1
    spawncenter = getrandompos(cellmatrix, cell_len)
    area, bb = initcellpositions!(cellmatrix, sigma, spawncenter, cell_len)
    edges = getedges(bb, cellmatrix)
    attrset = AttrSet(
        sigma,
        veg,
        area,
        getcenter(bb),
        bb
    )
    push!(cells.attrsets, attrset)
    append!(cells.edges, edges)
end

function setup(;ncells::Integer,
               size::Integer,
               cell_length::Integer)
    cells = Cells()
    cellmatrix = zeros(Int, size, size)
    for _ in 1:ncells
        spawncell!(cells, cellmatrix, cell_length)
    end
    cells, cellmatrix
end