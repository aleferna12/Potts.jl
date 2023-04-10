using GLMakie
GLMakie.activate!(fxaa=false)

function plotsimulation!(cellfield::Matrix, filename::AbstractString)
    scene = Scene(camera=campixel!, resolution=size(cellfield))
    image!(scene, cellfield, colormap=:nipy_spectral)
    save(filename, scene)
end
