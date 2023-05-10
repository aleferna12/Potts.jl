module CPM

export main

using Dates
using Images
using ImageView

using Statistics: mean
using Random: randperm, seed!

import Gtk

const PROJECT_HOME = dirname(dirname(@__FILE__))

include("utils.jl")
include("boundingbox.jl")
include("cell.jl")
include("parameters.jl")
include("environment.jl")
include("io.jl")

"Entry point."
function main(args=ARGS)
    params = length(args) > 0 ? readparams(Parameters, args[1], args[2:end]) : Parameters()
    runsimulation(Dish, Cells, params)
end

function runsimulation(envtype::Type{ET}, cellstype::Type{CT}, params) where {ET <: Environment, CT <: AbstractCells}
    makesimdirs(params.simdir, [params.imagesdirname], params.replaceprevsim)

    outputobjs = setupoutput(params)
    env = setupenv(envtype, cellstype, params)
    for i in 0:params.endsim
        step(env, i, params)
        if i % params.outputperiod == 0
            output!(env, i, outputobjs..., params)
    end end
end
 
end # module CPM
