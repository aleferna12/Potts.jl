module CPM

export main

using Dates
using Images
using ImageView

using Statistics: mean
using Random: randperm, seed!

import Gtk

include("utils.jl")
include("boundingbox.jl")
include("cell.jl")
include("parameters.jl")
include("dish.jl")
include("io.jl")

"Entry point."
function main(args=ARGS)
    params = length(args) == 1 ? readparams(args[1]) : Parameters()
    makesimdirs(params.simdir, [params.imagesdirname], params.replaceprevsim)
    
    outputobjs = setupoutput(params)
    cells = setup(params)
    for i in 0:params.endsim
        step(cells, i, params)
        if i % params.outputperiod == 0
            output!(cells, i, outputobjs..., params)
        end
end end
 
end # module CPM
