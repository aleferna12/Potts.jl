"A position on the simulation matrix."
struct Pos{T<:Number}
    x::T
    y::T
end
getx(pos::Pos) = pos.x
gety(pos::Pos) = pos.y

"Often in CPM simulations positions are represented as matrix indices, so we provide this alias definition."
const MatrixPos = Pos{Int}
Base.getindex(array::Array, pos::MatrixPos) = array[pos.x, pos.y]
Base.setindex!(array::Array, value, pos::MatrixPos) = array[pos.x, pos.y] = value
Base.checkindex(::Type{Bool}, inds::AbstractUnitRange, index::MatrixPos) = 0 < index.x * index.y <= length(inds)
getadjacentx(pos::MatrixPos) = [MatrixPos(pos.x - 1, pos.y), MatrixPos(pos.x + 1, pos.y)]
getadjacenty(pos::MatrixPos) = [MatrixPos(pos.x, pos.y - 1), MatrixPos(pos.x, pos.y + 1)]
vonneumann_neighbors(pos::MatrixPos) = [getadjacentx(pos); getadjacenty(pos)]
moore_neighbors(pos::MatrixPos) = [vonneumann_neighbors(pos); [MatrixPos(pos.x - 1, pos.y - 1),
                                                               MatrixPos(pos.x - 1, pos.y + 1),
                                                               MatrixPos(pos.x + 1, pos.y - 1),
                                                               MatrixPos(pos.x + 1, pos.y + 1)]]

function getrandompos(matrix::Matrix, borderpadding::Integer=0)
    xrange = range(1 + borderpadding, size(matrix)[1] - borderpadding)
    yrange = range(1 + borderpadding, size(matrix)[2] - borderpadding)
    MatrixPos(rand(xrange), rand(yrange))
end


abstract type AbstractTimer end
"Amount of time that must go by for a timer to fire."
getcooldown(timer::AbstractTimer) = timer.cooldown
"What time was the timer last reset (either automatically or manually)."
getlasttime(timer::AbstractTimer) = timer.lasttime
"How much time has elapsed between the last call to 'fire!()' or 'tick!()' 
and the last time the timer was reset (either automatically or manually)."
getelapsedtime(timer::AbstractTimer) = timer.elapsedtime
"Whether the timer resets automatically after firing."
autoresets(timer::AbstractTimer) = timer.autoreset
isactive(timer::AbstractTimer) = timer.active
setactive!(timer::AbstractTimer, val) = timer.active = val
reset!(timer::AbstractTimer, now) = begin timer.elapsedtime = 0; timer.lasttime = now end
"Sets the couter to fire next 'fire!()' or 'tick!()' call."
set!(timer::AbstractTimer) = timer.lasttime -= getcooldown(timer)

"Checks if it's time for the timer to fire.

Returns 'true' if the elapsed time between 'now' and the last reset of the timer has surpassed the timer's cooldown."
function fire!(timer::AbstractTimer, now)
    elapsed = now - getlasttime(timer)
    timer.elapsedtime = elapsed
    if isactive(timer) && elapsed >= getcooldown(timer)
        if autoresets(timer)
            reset!(timer, now)
        end
        return true
    end
    return false
end

"Updates the timer with the elapsed time since last update and checks if it's time to fire.

Returns 'true' if the total elapsed time since the last reset of the timer (including this update's 'elapsed') has surpassed the timer's cooldown."
function tick!(timer::AbstractTimer, elapsed)
    fire!(timer, getlasttime(timer) + getelapsedtime(timer) + elapsed)
end

mutable struct IterationTimer <: AbstractTimer
    cooldown::Int64
    elapsedtime::Int64
    lasttime::Int64
    autoreset::Bool
    active::Bool
end

function IterationTimer(cooldown; 
                 set=false, 
                 elapsedtime=0, 
                 lasttime=0, 
                 autoreset=true, 
                 active=true)
    timer = IterationTimer(cooldown, elapsedtime, lasttime, autoreset, active)
    if set
        set!(timer)
    end
    timer
end

"Handy method for 'tick!(timer, 1)' that is meant to be used every iteration of a loop."
tick!(timer::IterationTimer) = tick!(timer, 1)

mutable struct Timer <: AbstractTimer
    cooldown::Int64
    elapsedtime::Int64
    lasttime::Int64
    autoreset::Bool
    active::Bool
    Timer(cooldown, elapsedtime, lasttime, autoreset, active) = new(in_ns(cooldown), in_ns(elapsedtime), in_ns(lasttime), autoreset, active)
end

function Timer(cooldown;
               set=false, 
               elapsedtime=0, 
               lasttime=time_ns(), 
               autoreset=true, 
               active=true) 
    timer = Timer(cooldown, elapsedtime, lasttime, autoreset, active)
    if set
        set!(timer)
    end
    timer
end

reset!(timer::Timer) = reset!(timer, time_ns())
fire!(timer::Timer, now::Period) = fire!(timer, in_ns(now))
"Handy method for 'fire!(timer, time_ns())'."
fire!(timer::Timer) = fire!(timer, time_ns())
tick!(timer::Timer, elapsed::Period) = tick!(timer, in_ns(elapsed))

in_ns(time::Real) = time
in_ns(time) = Dates.tons(time)


const Edge = Pair{MatrixPos, MatrixPos}

function removeedges!(edgeset, pos1, pos2)
    pop!(edgeset, Edge(pos1, pos2), Edge(pos2, pos1))
end

function addedges!(edgeset, pos1, pos2)
    push!(edgeset, Edge(pos1, pos2), Edge(pos2, pos1))
end


allsame(x) = all(y -> y == first(x), x)

ordered(x1, x2) = x2 < x1 ? (x2, x1) : (x1, x2)

function fullyconcrete(type::Type)
    stack = Type[type]
    while !isempty(stack)
        t = pop!(stack)
        if !isconcretetype(t)
            return false
        end
        paramtypes = filter(pt -> pt isa Type, collect(t.parameters))
        append!(stack, fieldtypes(t), paramtypes)
    end
    true
end