module CPM

export main

using Images

using Statistics: mean
using Random: randperm, seed!

include("utils.jl")
include("boundingbox.jl")
include("cell.jl")
include("parameters.jl")
include("dish.jl")
include("io.jl")

"Entry point."
function main()
    cells = setup()
    plotsimulation!(cells, "run/t0.png")
    for i in 1:PARAMS.endsim
        step(i, cells)
        if i % PARAMS.imageperiod == 0
            println("Timestep: $i")
            plotsimulation!(cells, "run/t$i.png")
    end end
end
 
end # module CPM
