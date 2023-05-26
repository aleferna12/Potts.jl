module CPM

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
include("parameters.jl")
include("environment.jl")
include("io.jl")

function Base.run(env::Environment, params, typechecks=true)
    if typechecks
        if !fullyconcrete(typeof(env))
            @warn "'$env' may contain abstract fields in it's type hierarchy"
        end
        if !fullyconcrete(typeof(params))
            @warn "'$params' may contain abstract fields in it's type hierarchy"
        end
    end

    starttime = now()
    makesimdirs(params.simdir, [params.imagesdirname], params.replaceprevsim)

    setupenv!(env, params)
    outputobjs = setupoutput(params)
    for i in 0:params.endsim
        step!(env, i, params)
        if i % params.outputperiod == 0
            output!(env, i, outputobjs..., params)
    end end
    endmessage(env, starttime)
end
 
end # module CPM
