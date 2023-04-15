@enum Tau veg div

mutable struct AttrSet
    sigma::Int
    tau::Tau
    area::Int
    center::Pos
    bb::BoundingBox
end

mutable struct Genome{}
    sigma::Int
end

struct Cells
    attrsets::Vector{AttrSet}
    genomes::Vector{Genome}
    edges::Vector{Edge}
end
Cells() = Cells(AttrSet[], Genome[], Edge[])
Base.length(cells::Cells) = length(cells.attrsets)
Base.lastindex(cells::Cells) = lastindex(cells.attrsets)