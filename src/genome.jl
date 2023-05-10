struct Receptor
    value::Float64
    scale::Float64
end
read(receptor::Receptor) = receptor.value
write!(receptor::Receptor, val) = receptor.value = val * receptor.scale

struct Gene
    active::Bool
    threshold::Float64
end

struct Genome
    receptors::Vector{Receptor}
    reggenes::Vector{Gene}
    outgenes::Vector{Gene}
    weights::Vector{Float64}
end
# TODO
# Genome(nin::Integer, nreg::Integer, nout::Integer) = ...