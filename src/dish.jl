using Statistics: mean

function randompos(field::Matrix, borderpadding=0)
    irange = range(1 + borderpadding, size(field)[1] - borderpadding)
    jrange = range(1 + borderpadding, size(field)[2] - borderpadding)
    rand(irange), rand(jrange)
end

function initcellpositions!(sigma::Int, i::Int, j::Int, cell_len::Int, cellfield::Matrix)
    area = 0
    bb = BoundingBox(i, i, j, j)
    rangebb = BoundingBox(
        i - floor(cell_len / 2), 
        i + ceil(cell_len / 2),
        j - floor(cell_len / 2),
        j + ceil(cell_len / 2)
    )
    for (i, j) in iteratepos(rangebb)
        if cellfield[i, j] == -1
            cellfield[i, j] = sigma
            area += 1
            bb = addpos!(bb, i, j)
        end
    end
    area, bb
end

"""Get the center of a cell."""
function cellcenterpos(cellfield::Matrix, sigma::Int, cellbb::BoundingBox)
    is, js = Int[], Int[]
    for (i, j) in iteratepos(cellbb)
        if cellfield[i, j] == sigma
            push!(is, i)
            push!(js, j)
        end
    end
    mean(is), mean(js)
end

function spawncell!(cells::Cells, cellfield::Matrix, cell_len::Int=5)
    sigma = length(cells)
    i, j = randompos(cellfield, cell_len + 3)  # +3 for safety
    area, cellbb = initcellpositions!(sigma, i, j, cell_len, cellfield)
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

function start_simulation!(ncells::Int)   
    for _ in 1:ncells
        spawncell!(cells, cellfield)
    end
end

cellfield = zeros(100, 100)
cells = Cells()