module CPM

export main

using Images

using Statistics: mean
using Random: randperm

include("utils.jl")
include("boundingbox.jl")
include("cell.jl")
include("dish.jl")
include("io.jl")

"Entry point."
function main()
    params = (
        ncells = 50,
        size = 100,
        cell_length = 5,
        boltzmanntemp = 16,
        targetarea = 100,
        sizelambda = 1,
        adhesiontable = zeros(2, 2),
        adhesionmedium = 0
    )
    cells = setup(;params...)
    for i in 1:10000
        step(i, cells; params...)
        if i % 1000 == 0
            plotsimulation!(cells, "run/t$i.png")
    end end
end
 
end # module CPM
