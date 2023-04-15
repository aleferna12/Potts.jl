using GLMakie
GLMakie.activate!(fxaa=false)

function plotsimulation!(cellmatrix::Matrix, filename::AbstractString)
    scene = Scene(camera=campixel!, resolution=size(cellmatrix))
    image!(scene, cellmatrix, colormap=:nipy_spectral)
    save(filename, scene)
end
