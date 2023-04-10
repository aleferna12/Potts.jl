include("CPM.jl")
import .CPM

cells, cellfield = CPM.setup(ncells=20,
                             size=100,
                             cell_length=3)
CPM.plotsimulation!(cellfield, "run/t0.png")