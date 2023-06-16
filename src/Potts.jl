module Potts

export run, defaultparams, readparams, CPM, Dish, Cell, EvolvableCell

using Dates
using Images
using ImageView

using Statistics: mean
using Random: randperm, seed!

import Gtk

const PROJECT_HOME = dirname(dirname(@__FILE__))

include("constants.jl")
include("utils.jl")
include("timer.jl")
include("pos.jl")
include("boundingbox.jl")
include("genome.jl")
include("cell.jl")
include("environment.jl")
include("model.jl")
include("io.jl")
include("parameters.jl")
include("examples/evolution.jl")

function Base.run(model::AbstractCPM)
    seed!(model[:seed])
    starttime = now()
    makesimdirs(model[:simdir], [model[:imagesdirname]], model[:replaceprevsim])

    setup!(model)
    outputobjs = setupoutput(;getparams(model)...)
    for i in 0:model[:endsim]
        step!(model)
        output!(model, i, outputobjs...)
    end
    endmessage(model, starttime)
end
 
end # module Potts
