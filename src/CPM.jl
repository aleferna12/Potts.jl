module CPM

export main

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
function main()
    makesimdirs(PARAMS.simdir, [PARAMS.imagesdirname], PARAMS.replaceprevsim)
    
    cells = setup()
    gui = makegui(
        PARAMS.displayperiod,
        length(PARAMS.imageplots),
        PARAMS.displaysize
    )

    for i in 0:PARAMS.endsim
        step(cells, i)
        output(cells, gui, i)
end end
 
end # module CPM
