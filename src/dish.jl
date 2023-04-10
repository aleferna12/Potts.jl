using Statistics: mean

function randompos(field::Matrix, borderpadding::Integer=0)
    irange = range(1 + borderpadding, size(field)[1] - borderpadding)
    jrange = range(1 + borderpadding, size(field)[2] - borderpadding)
    rand(irange), rand(jrange)
end

function initcellpositions!(cellfield::Matrix, sigma::Integer, i::Integer, j::Integer, cell_len::Integer)
    area = 0
    bb = BoundingBox(i, i, j, j)
    rangebb = BoundingBox(
        i - floor(cell_len / 2), 
        i + ceil(cell_len / 2) - 1,
        j - floor(cell_len / 2),
        j + ceil(cell_len / 2) - 1
    )
    for (i, j) in iteratepos(rangebb)
        if cellfield[i, j] == 0
            cellfield[i, j] = sigma
            area += 1
            addpos!(bb, i, j)
        end
    end
    area, bb
end

"""Get the center of a cell."""
function cellcenterpos(cellfield::Matrix, sigma::Integer, cellbb::BoundingBox)
    is, js = Int[], Int[]
    for (i, j) in iteratepos(cellbb)
        if cellfield[i, j] == sigma
            push!(is, i)
            push!(js, j)
        end
    end
    mean(is), mean(js)
end

function spawncell!(cells::Cells, cellfield::Matrix, cell_len::Integer)
    sigma = lastindex(cells) + 1
    i, j = randompos(cellfield, cell_len)
    area, cellbb = initcellpositions!(cellfield, sigma, i, j, cell_len)
    x, y = cellcenterpos(cellfield, sigma, cellbb)
    attrset = AttrSet(
        sigma,
        veg,
        area,
        x, y,
        cellbb
    )
    push!(cells.attrsets, attrset)
end

function setup(;ncells::Integer,
               size::Integer,
               cell_length::Integer)
    cells = Cells()
    cellfield = zeros(size, size)
    for _ in 1:ncells
        spawncell!(cells, cellfield, cell_length)
    end
    cells, cellfield
end