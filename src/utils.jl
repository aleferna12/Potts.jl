"A position on the simulation matrix."
struct Pos{T<:Number}
    x::T
    y::T
end
getx(pos::Pos) = pos.x
gety(pos::Pos) = pos.y
Base.round(postype::Type, pos::Pos) = postype(round(getx(pos)), round(gety(pos)))
Base.round(pos::Pos) = round(typeof(pos), pos)
Base.convert(postype::Type{<:Pos}, pos::NTuple{2}) = postype(pos[1], pos[2])

"Determines on which side of a line (specified by 'm' and 'n') a point lies. Returns 1 for 'up' or 'left' and -1 otherwise."
function whichside(m, n, pos::Pos)
    y = m * getx(pos) + n
    sign(gety(pos) - y)
end

"Often in CPM simulations positions are represented as matrix indices, so we provide this alias definition."
const MatrixPos = Pos{Int}
Base.getindex(matrix::Matrix, pos::MatrixPos) = matrix[pos.x, pos.y]
Base.getindex(matrix::Matrix, positions::MatrixPos...) = [matrix[pos] for pos in positions]
Base.setindex!(matrix::Matrix, value, pos::MatrixPos) = matrix[pos.x, pos.y] = value
Base.setindex!(matrix::Matrix, value, positions::MatrixPos...) = for pos in positions matrix[pos] = value end
Base.checkbounds(::Type{Bool}, matrix::Matrix, index::MatrixPos) = checkbounds(Bool, matrix, CartesianIndex(getx(index), gety(index)))
Base.convert(::Type{MatrixPos}, pos::Pos) = MatrixPos(getx(pos), gety(pos))
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
rewind!(timer::IterationTimer, val) = timer.elapsed -= val

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