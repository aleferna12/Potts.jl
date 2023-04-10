import Base: length, lastindex

@enum Tau veg div

mutable struct AttrSet
    sigma::Int
    tau::Tau
    area::Int
    x::Float64
    y::Float64
    boundingbox::BoundingBox
end

mutable struct Genome{}
    sigma::Int
end

struct Cells
    attrsets::Vector{AttrSet}
    genomes::Vector{Genome}
end
Cells() = Cells(AttrSet[], Genome[])
length(cells::Cells) = length(cells.attrsets)
lastindex(cells::Cells) = lastindex(cells.attrsets)