"A position on the simulation matrix."
struct Pos{T<:Number}
    x::T
    y::T
end
getx(pos::Pos) = pos.x
gety(pos::Pos) = pos.y
inbounds(pos::Pos, bounds::Tuple) = 0 < pos.x <= bounds[1] && 0 < pos.y <= bounds[2]
filterinbounds(positions, bounds::Tuple) = filter(pos -> inbounds(pos, bounds), positions)

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


abstract type AbstractCounter end
"Amount of time that must go by for a counter to fire."
getcooldown(counter::AbstractCounter) = counter.cooldown
"What time was the counter last reset (either automatically or manually)."
getlasttime(counter::AbstractCounter) = counter.lasttime
"How much time has elapsed between the last call to 'fire!()' or 'tick!()' 
and the last time the counter was reset (either automatically or manually)."
getelapsedtime(counter::AbstractCounter) = counter.elapsedtime
"Whether the counter resets automatically after firing."
autoresets(counter::AbstractCounter) = counter.autoreset
isactive(counter::AbstractCounter) = counter.active
setactive!(counter::AbstractCounter, val) = counter.active = val
reset!(counter::AbstractCounter, now) = begin counter.elapsedtime = 0; counter.lasttime = now end
"Sets the couter to fire next 'fire!()' or 'tick!()' call."
set!(counter::AbstractCounter) = counter.lasttime -= getcooldown(counter)

"Checks if it's time for the counter to fire.

Returns 'true' if the elapsed time between 'now' and the last reset of the counter has surpassed the counter's cooldown."
function fire!(counter::AbstractCounter, now)
    elapsed = now - getlasttime(counter)
    counter.elapsedtime = elapsed
    if isactive(counter) && elapsed >= getcooldown(counter)
        if autoresets(counter)
            reset!(counter, now)
        end
        return true
    end
    return false
end

"Updates the counter with the elapsed time since last update and checks if it's time to fire.

Returns 'true' if the total elapsed time since the last reset of the counter (including this update's 'elapsed') has surpassed the counter's cooldown."
function tick!(counter::AbstractCounter, elapsed)
    fire!(counter, getlasttime(counter) + getelapsedtime(counter) + elapsed)
end

mutable struct Counter <: AbstractCounter
    cooldown::Int64
    elapsedtime::Int64
    lasttime::Int64
    autoreset::Bool
    active::Bool
end

function Counter(cooldown; 
                 set=false, 
                 elapsedtime=0, 
                 lasttime=0, 
                 autoreset=true, 
                 active=true)
    counter = Counter(cooldown, elapsedtime, lasttime, autoreset, active)
    if set
        set!(counter)
    end
    counter
end

"Handy method for 'tick!(counter, 1)' that is meant to be used every iteration of a loop."
tick!(counter::Counter) = tick!(counter, 1)

mutable struct Timer <: AbstractCounter
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

orderedpair(x1, x2) = x2 < x1 ? x2 => x1 : x1 => x2