abstract type AbstractCell end
getsigma(cell::AbstractCell) = cell.sigma
gettau(cell::AbstractCell) = cell.tau
getcenter(cell::AbstractCell) = cell.center
setcenter!(cell::AbstractCell, val) = cell.center = val
getarea(cell::AbstractCell) = cell.area
setarea!(cell::AbstractCell, val) = cell.area = val
addarea!(cell::AbstractCell, val) = setarea!(cell, getarea(cell) + val)


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

mutable struct Cell <: AbstractCell
    sigma::Int
    tau::Int
    area::Int
    center::Pos{Float64}
end