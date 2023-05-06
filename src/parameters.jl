"Parameters of the simulation that can be accessed in the global scope by highly parametric functions."
Base.@kwdef struct Parameters
    endsim::Int = 1000000
    fieldsize::Int = 100
    ncells::Int = 15
    cell_length::Int = 5
    boltzmanntemp::Float64 = 12
    sizelambda::Float64 = 1
    targetcellarea::Int = 25  # TODO: change to float (also on function calls)
    adhesiontable::Matrix{Float64} = fill(8, 2, 2)
    adhesionmedium::Float64 = 8
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
end

function readparams(filepath)
    content = read(filepath, String)
    parsable = replace(content, r"::.*(?==)" => "") # The type annotations on the file are decorative what matters is the struct types
    parsed = Meta.parse("begin\n" * parsable * "\nend")
    paramtup = NamedTuple(expr.args[1] => eval(expr.args[2]) for expr in parsed.args if expr isa Expr)
    Parameters(; paramtup...)
end