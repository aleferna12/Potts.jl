mutable struct CellAttrs
    sigma::Int
    tau::Int
    area::Int
    center::Pos{Float64}
end

abstract type AbstractCells end
# This interface is designed to automatically work with subtypes of AbstractCells, provided that they implement a cellattrs field
# A subtype of AbstractCells that doesn't implement this field will still work as long as it instead implements new accessor methods
getsigmas(cells::AbstractCells) = [cells.cellattrs[i].sigma for i in 1:length(cells)]
gettau(cells::AbstractCells, sigma) = cells.cellattrs[sigma].tau
gettaus(cells::AbstractCells) = [cells.cellattrs[i].tau for i in 1:length(cells)]
getcenter(cells::AbstractCells, sigma) = cells.cellattrs[sigma].center
getcenters(cells::AbstractCells) = [cells.cellattrs[i].center for i in 1:length(cells)]
getarea(cells::AbstractCells, sigma) = cells.cellattrs[sigma].area
getareas(cells::AbstractCells) = [cells.cellattrs[i].area for i in 1:length(cells)]
addarea!(cells::AbstractCells, sigma, val) = cells.cellattrs[sigma].area += val
Base.length(cells::AbstractCells) = length(cells.cellattrs)
Base.lastindex(cells::AbstractCells) = lastindex(cells.cellattrs)
Base.size(cells::AbstractCells) = size(getmatrix(cells))
Base.size(cells::AbstractCells, dim::Integer) = size(getmatrix(cells), dim)

"Holds information about a group of cells. Uses the components pattern to miniminze cache misses and improve performance,
although this implementation has a single component (cellattrs)."
struct Cells <: AbstractCells
    cellattrs::Vector{CellAttrs}
end
Cells() = Cells(CellAttrs[])
Base.push!(cells::Cells, sigma, tau, area, center) = push!(cells.cellattrs, CellAttrs(sigma, tau, area, center))

struct EvolvableCells <: AbstractCells
    cellattrs::Vector{CellAttrs}
    genomes::Vector{Genome}
end

function getadhenergy(cells::AbstractCells, sigma1, sigma2, adhesiontable, adhesionmedium, borderenergy)::Float64
    if sigma1 == sigma2
        return 0.
    end
    sigma1, sigma2 = ordered(sigma1, sigma2)
    if sigma1 < 0
        return sigma2 == 0 ? 0. : borderenergy
    elseif sigma1 == 0
        return adhesionmedium[gettau(cells, sigma2)]
    end
    adhesiontable[gettau(cells, sigma1), gettau(cells, sigma2)] * 2
end