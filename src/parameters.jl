abstract type AbstractParameters end

"Parameters of the simulation."
Base.@kwdef struct Parameters <: AbstractParameters
    endsim::Int = 1000000
    fieldsize::Int = 100
    taus::Int = 1
    ncells::Vector{Int} = [20]
    cell_length::Vector{Int} = [5]
    targetcellarea::Vector{Float64} = [50]
    adhesiontable::Matrix{Float64} = hcat(8)
    adhesionmedium::Vector{Float64} = [8]
    cellcolors::Vector{RGB} = [colorant"gray"]
    boltzmanntemp::Float64 = 12
    sizelambda::Float64 = 1
    borderenergy::Float64 = 100
    simdir::String = "./runs/debug"
    imagesdirname::String = "image"
    replaceprevsim::Bool = false
    imageplots::Vector{Symbol} = [:sigma]
    savegif::Bool = false
    outputperiod::Int = 10
    infoperiod::Int = 1000
    imageperiod::Int = 1000
    displayframerate::Int = 30
    displaysize::Int = 350
    drawcellborders::Bool = true
    drawcellcenters::Bool = false
end

"Read a parameters file into a type 'paramstype'. This type must implement a keyword constructor (using the 'Base.@kwdef' macro)."
function readparams(paramstype::Type{T}, filepath, args) where T <: AbstractParameters
    content = readresource(filepath) * '\n' * join(args, '\n')
    parsable = replace(content, r"::.*(?==)" => "") # The type annotations on the file are decorative what matters is the struct types
    parsed = Meta.parse("begin\n" * parsable * "\nend")
    paramtup = NamedTuple(expr.args[1] => eval(expr.args[2]) for expr in parsed.args if expr isa Expr)
    paramstype(; paramtup...)
end

function readresource(filepath)
    try
        return read(filepath, String)
    catch e
        respath = joinpath(PROJECT_HOME, "resources")
        if e isa SystemError && filepath in readdir(respath)
            # File not found, try to find it in resources
            return read(joinpath(respath, filepath), String)
        else
            rethrow()
end end end