"Parameters of the simulation that can be accessed in the global scope by highly parametric functions."
Base.@kwdef struct Parameters0
    "Size of the field."
    fieldsize::Int = 100
    "Number of starting cells."
    ncells::Int = 15
    "Starting cell length."
    cell_length::Int = 5
    "Target area for cells."
    targetcellarea::Int = 25
    "Size constraint scale (how strongly cells deviating from 'targetcellarea' are penalized)."
    sizelambda::Float64 = 1
    "Symmetric table used to determine cell-cell adhesion energy. The indices are the taus of cells."
    adhesiontable::Matrix{Float64} = fill(16, 2, 2)
    "Adhesion energy of a cell with the medium."
    adhesionmedium::Float64 = 8
    "Energy on the interface with a border site."
    borderenergy::Float64 = 100
    "Temperature controling the 'thermodynamic' behaviour of the system."
    boltzmanntemp::Float64 = 16
    "How many MCSs to run the simulation for."
    endsim::Int = 1000000
    "Directory where the simulation files will be created."
    simdir::String = "./run"
    "Name of the subdirectory where image files are stored."
    imagesdirname::String = "image"
    "Whether to delete the directory indicated by 'simdir' if it already exists."
    replaceprevsim::Bool = false
    "Vector of which real-time plots to make of the simulation. Possible values are: sigma, tau."
    imageplots::Vector{Symbol} = [:sigma]
    "Whether to save GIFs of the simulation instead of individual image files."
    savegif::Bool = true
    "How often in MCSs to display information about the simulation on the console. Non-positive values disable this feature."
    infoperiod::Int = 1000
    "How often in MCSs we save an image (or frame if 'savegif' is true) of the simulation. Non-positive values disable this feature."
    imageperiod::Int = 1000
    "How often in MCSs the real-time display of the simulation should be updated. Non-positive values disable this feature."
    displayperiod::Int = 100
    "How big is the window of the real-time display."
    displaysize::Int = 350
    "Whether to draw borders between cells."
    drawcellborders::Bool = true
end

Parameters = Parameters0 # This is a trick to facilitate use of Revise, see docs
# This should be a const but Revise doesnt like consts so im keeping it like that for now 
PARAMS = Parameters(  # TODO: Replace with readparameters function
    endsim=10000,
    imageperiod=-1,
    # displayperiod=-1,
    imageplots=[:sigma, :tau],
    displaysize=350,
    replaceprevsim=true,
)