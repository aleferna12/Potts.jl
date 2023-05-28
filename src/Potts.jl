module Potts

export run

using Dates
using Images
using ImageView

using Statistics: mean
using Random: randperm, seed!

import Gtk

const PROJECT_HOME = dirname(dirname(@__FILE__))

include("utils.jl")
include("boundingbox.jl")
include("genome.jl")
include("cell.jl")
include("environment.jl")
include("model.jl")
include("io.jl")
include("parameters.jl")

function Base.run(model::AbstractCPM)
    seed!(model[:seed])
    starttime = now()
    makesimdirs(model[:simdir], [model[:imagesdirname]], model[:replaceprevsim])

    setup!(model)
    outputobjs = setupoutput(;getparams(model)...)
    for i in 0:model[:endsim]
        step!(model, i)
        if i % model[:outputperiod] == 0
            output!(model, i, outputobjs...)
    end end
    endmessage(model, starttime)
end
 
end # module Potts
