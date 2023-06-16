"Abstract type representing a cell. This interface assumes some fields are implemented."
abstract type AbstractCell end
# TODO: try getting rid of the abstract getters and writing a macro similar to Lazy.@forward
#       then measure performance against current implementation
getsigma(cell::AbstractCell) = cell.sigma
gettau(cell::AbstractCell) = cell.tau
getcenter(cell::AbstractCell) = cell.center
setcenter!(cell::AbstractCell, val) = cell.center = val
getarea(cell::AbstractCell) = cell.area
setarea!(cell::AbstractCell, val) = cell.area = val
addarea!(cell::AbstractCell, val) = setarea!(cell, getarea(cell) + val)
gettargetarea(cell::AbstractCell) = cell.targetarea
settargetarea!(cell::AbstractCell, val) = cell.targetarea = val
isalive(cell::AbstractCell) = cell.alive
kill!(cell::AbstractCell) = cell.alive = false
growthperc(cell::AbstractCell) = getcooldown(cell.divtimer) == 0 ? 0. : getelapsed(cell.divtimer) / getcooldown(cell.divtimer)
isdividing(cell::AbstractCell) = isactive(cell.divtimer)
getdivtime(cell::AbstractCell) = getcooldown(cell.divtimer)
dividenow!(cell::AbstractCell) = fire!(cell.divtimer)
startdividing!(cell::AbstractCell) = begin activate!(cell.divtimer); reset!(cell.divtimer) end
stopdividing!(cell::AbstractCell) = deactivate!(cell.divtimer)

function addmomentum!(cell::AbstractCell, vec::Pos)
    center = getcenter(cell)
    newcenter = Pos(getx(center) + (getx(vec) - getx(center)) / getarea(cell),
                    gety(center) + (gety(vec) - gety(center)) / getarea(cell))
    setcenter!(cell, newcenter)
end

function removemomentum!(cell::AbstractCell, vec::Pos)
    center = getcenter(cell)
    newcenter = Pos(getx(center) - (getx(vec) - getx(center)) / getarea(cell),
                    gety(center) - (gety(vec) - gety(center)) / getarea(cell))
    setcenter!(cell, newcenter)
end

function gainpos!(cell::AbstractCell, pos::MatrixPos)
    addarea!(cell, 1)
    addmomentum!(cell, pos)
end

function losepos!(cell::AbstractCell, pos::MatrixPos)
    addarea!(cell, -1)
    removemomentum!(cell, pos)
end

abstract type EvolvableCell <: AbstractCell end

Base.@kwdef mutable struct Cell <: AbstractCell
    sigma::Int
    tau::Int
    area::Int
    targetarea::Float64
    center::Pos{Float64}
    alive::Bool = true
    divtimer::IterationTimer
end