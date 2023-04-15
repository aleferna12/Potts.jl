using Statistics: mean

function setup(;ncells::Integer,
               size::Integer,
               cell_length::Integer)
    cells = Cells(size)
    for _ in 1:ncells
        spawncell!(cells, cell_length)
    end
    cells
end