include("CPM.jl")
import .CPM

cells, cellmatrix = CPM.setup(ncells=20,
                              size=100,
                              cell_length=5)
CPM.plotsimulation!(cellmatrix, "run/t0.png")
