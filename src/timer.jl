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
    false
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