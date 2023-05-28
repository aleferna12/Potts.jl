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
moore_neighbors(pos::MatrixPos) = [MatrixPos(pos.x - 1, pos.y),
                                   MatrixPos(pos.x + 1, pos.y),
                                   MatrixPos(pos.x, pos.y - 1),
                                   MatrixPos(pos.x, pos.y + 1),
                                   MatrixPos(pos.x - 1, pos.y - 1),
                                   MatrixPos(pos.x - 1, pos.y + 1),
                                   MatrixPos(pos.x + 1, pos.y - 1),
                                   MatrixPos(pos.x + 1, pos.y + 1)]

function getrandompos(matrix::Matrix, borderpadding::Integer=0)
    xrange = range(1 + borderpadding, size(matrix)[1] - borderpadding)
    yrange = range(1 + borderpadding, size(matrix)[2] - borderpadding)
    MatrixPos(rand(xrange), rand(yrange))
end


abstract type AbstractTimer end
"Whether the timer resets automatically after firing."
autoresets(timer::AbstractTimer) = timer.autoreset
isactive(timer::AbstractTimer) = timer.active
activate!(timer::AbstractTimer) = timer.active = true
deactivate!(timer::AbstractTimer) = timer.active = false

mutable struct IterationTimer <: AbstractTimer
    cooldown::Int64
    elapsed::Int64
    autoreset::Bool
    active::Bool
end
function IterationTimer(cooldown; 
                        set=false, 
                        elapsedtime=0,
                        autoreset=true, 
                        active=true)
    timer = IterationTimer(cooldown, elapsedtime, autoreset, active)
    if set
        set!(timer)
    end
    timer
end
"Amount of time that must go by for a timer to fire."
getcooldown(timer::IterationTimer) = timer.cooldown
"How many iterations have elapsed between the last call to 'fire!()' 
and the last time the timer was reset (either automatically or manually)."
getelapsed(timer::IterationTimer) = timer.elapsed
reset!(timer::IterationTimer) = timer.elapsed = 0
"Sets the timer to fire next 'fire!()' call."
set!(timer::IterationTimer) = timer.elapsed = getcooldown(timer)

"Returns 'true' if it's time for the timer to fire and 'false' otherwise."
function fire!(timer::IterationTimer, elapsed=1)
    timer.elapsed += elapsed
    if isactive(timer) && getelapsed(timer) >= getcooldown(timer)
        if autoresets(timer)
            reset!(timer)
        end
        return true
    end
    return false
end

mutable struct Timer <: AbstractTimer
    cooldown::Int64
    lastreset::Int64
    autoreset::Bool
    active::Bool
    Timer(cooldown, lastreset, autoreset, active) = new(in_ns(cooldown), in_ns(lastreset), autoreset, active)
end
function Timer(cooldown;
               set=false, 
               lastreset=time_ns(), 
               autoreset=true, 
               active=true) 
    timer = Timer(cooldown, lastreset, autoreset, active)
    if set
        set!(timer)
    end
    timer
end
"Amount of time that must go by for a timer to fire."
getcooldown(timer::Timer) = timer.cooldown
"What time was the timer last reset (either automatically or manually)."
getlastreset(timer::Timer) = timer.lastreset

reset!(timer::Timer, now=time_ns()) = timer.lastreset = now
"Sets the timer to fire next 'fire!()' call."
set!(timer::Timer, now=time_ns()) = timer.lastreset = now - getcooldown(timer)

"Returns 'true' if it's time for the timer to fire and 'false' otherwise."
function fire!(timer::Timer, now=time_ns())
    if isactive(timer) && now - getlastreset(timer) >= getcooldown(timer)
        if autoresets(timer)
            reset!(timer, now)
        end
        return true
    end
    return false
end

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