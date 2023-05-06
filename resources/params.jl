###########
# General #
###########

# How many MCSs to run the simulation for
endsim::Int = 1000000

# Size of the field
fieldsize::Int = 100

#######################
# Starting Conditions #
#######################

# Number of starting cells
ncells::Int = 15

# Starting cell length
cell_length::Int = 5

###############
# Hamiltonian #
###############

# Temperature controling the 'thermodynamic' behaviour of the system
boltzmanntemp::Float64 = 12

# Size constraint scale (how strongly cells deviating from 'targetcellarea' are penalized)
sizelambda::Float64 = 1

# Target area for cells
targetcellarea::Int = 25

# Symmetric table used to determine cell-cell adhesion energy
# The indices are the taus of cells
adhesiontable::Matrix{Float64} = fill(8, 2, 2)

# Adhesion energy of a cell with the medium
adhesionmedium::Float64 = 8

# Energy on the interface with a border site
borderenergy::Float64 = 100

##########
# Output #
##########

# Directory where the simulation files will be created
simdir::String = "./run"

# Whether to delete the directory indicated by 'simdir' if it already exists
replaceprevsim::Bool = false

# Name of the subdirectory where image files are stored
imagesdirname::String = "image"

# Vector of which real-time plots to make of the simulation
# Possible values are: sigma, tau
imageplots::Vector{Symbol} = [:sigma]

# Whether to save GIFs of the simulation instead of individual image files
# Not recommended for long simulations (memory might fill up)
savegif::Bool = false

# How many MCSs must have passed before we try to save any output of the simulation
# Can be used to remove overhead of output functions
# Non-positive values disable this feature
outputperiod::Int = 10

# How often in MCSs (+ 'outputperiod') to display information about the simulation on the console
# Non-positive values disable this feature
infoperiod::Int = 1000

# How often in MCSs, (+ 'outputperiod') we save an image (or frame if 'savegif' is true) of the simulation
# Non-positive values disable this feature
imageperiod::Int = 1000

# Framerate for the display updates (in real life time)
# Non-positive values disable this feature
# Low framerate may result in an unresponsive window
# This parameter is also affected by 'outputperiod'
displayframerate::Int = 30

# How big is the window of the real-time display
displaysize::Int = 350

# Whether to draw borders between cells
drawcellborders::Bool = true
