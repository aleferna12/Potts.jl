include("CPM.jl")
import .CPM

CPM.start_simulation!(100)
CPM.plotsimulation!(CPM.cellfield, "run/t0.png")