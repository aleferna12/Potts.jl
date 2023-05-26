"Read a parameters file into a NamedTuple."
function readparams(filepath)
    contents = readresource(filepath)
    parsable = Meta.parse("begin\n" * contents * "\nend")
    exprs = filter(expr -> !isa(expr, LineNumberNode), parsable.args)
    typednamedtuple(exprs)
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

macro typednamedtuple(ntexpr)
    :($(typednamedtuple(ntexpr.args)))
end

function typednamedtuple(exprs)
    pairs = []
    for expr in exprs
        symbol, type, val = assignmentelements(expr)
        push!(pairs, symbol => convert(type,  val))
    end
    NamedTuple(pairs)
end

function assignmentelements(expr)
    if expr.head !== :(=) || length(expr.args) != 2
        throw("'expr' must be an assignment expression")
    end

    assigned, value = expr.args[1], eval(expr.args[2])
    if assigned isa Symbol
        return assigned, typeof(value), value
    end
    symbol, type = assigned.args
    symbol, eval(type), value
end

"Default parameters of the simulation."
const DEFAULTPARAMETERS = @typednamedtuple (
    endsim::Int = 100000,
    fieldsize::Int = 100,
    taus::Int = 1,
    ncells::Vector{Int} = [20],
    cell_length::Vector{Int} = [5],
    targetcellarea::Vector{Float64} = [50],
    adhesiontable::Matrix{Float64} = hcat(8),
    adhesionmedium::Vector{Float64} = [8],
    cellcolors::Vector{RGB{N0f8}} = [colorant"gray"],
    boltzmanntemp::Float64 = 12,
    sizelambda::Float64 = 1,
    borderenergy::Float64 = 100,
    simdir::String = "./runs/debug",
    imagesdirname::String = "image",
    replaceprevsim::Bool = false,
    imageplots::Vector{Symbol} = [:sigma],
    savegif::Bool = false,
    outputperiod::Int = 10,
    infoperiod::Int = 1000,
    imageperiod::Int = 1000,
    displayframerate::Int = 30,
    displaysize::Int = 350,
    drawcellborders::Bool = true,
    drawcellcenters::Bool = false
)