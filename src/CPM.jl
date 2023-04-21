module CPM

export main

using Images

using Statistics: mean
using Random: randperm

include("utils.jl")
include("parameters.jl")
include("boundingbox.jl")
include("cell.jl")
include("dish.jl")
include("io.jl")

"Entry point."
function main()
    cells = setup()
    for i in 1:10000
        step(i, cells)
        if i % 250 == 0
            println("Timestep: $i")
            plotsimulation!(cells, "run/t$i.png")
    end end
end
 
end # module CPM
