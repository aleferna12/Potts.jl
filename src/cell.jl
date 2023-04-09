@enum Tau mig div

mutable struct AttrSet
    sigma::Int
    tau::Tau
    area::Int
end

mutable struct Genome
    sigma::Int
end

struct Cells
    attribute_sets::Vector{AttrSet}
    genomes::Vector{Genome}
end
Cells() = Cells(AttrSet[], Genome[])