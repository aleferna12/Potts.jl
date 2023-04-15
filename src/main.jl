include("CPM.jl")
import .CPM

cells = CPM.setup(ncells=20,
                  size=100,
                  cell_length=5)
for e in cells.edgeset
    cells.matrix[e.pos1.x, e.pos1.y] -= 1
    cells.matrix[e.pos2.x, e.pos2.y] -= 1
end
CPM.plotsimulation!(cells.matrix, "run/t0.png")