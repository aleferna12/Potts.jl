###########
# General #
###########

# How many MCSs to run the simulation for
endsim::Integer = 1000000

# Size of the field
fieldsize::Integer = 100

#######################
# Starting Conditions #
#######################

# Number of starting cells
ncells::Integer = 15

# Starting cell length
cell_length::Integer = 5

###############
# Hamiltonian #
###############

# Temperature controling the 'thermodynamic' behaviour of the system
boltzmanntemp::Real = 12

# Size constraint scale (how strongly cells deviating from 'targetcellarea' are penalized)
sizelambda::Real = 1

# Target area for cells
targetcellarea::Integer = 25

# Symmetric table used to determine cell-cell adhesion energy
# The indices are the taus of cells
adhesiontable::Matrix{Real} = fill(8, 2, 2)

# Adhesion energy of a cell with the medium
adhesionmedium::Real = 8

# Energy on the interface with a border site
borderenergy::Real = 100

##########
# Output #
##########

# Directory where the simulation files will be created
simdir = "./run"

# Whether to delete the directory indicated by 'simdir' if it already exists
replaceprevsim = false

# Name of the subdirectory where image files are stored
imagesdirname = "image"

# Vector of which real-time plots to make of the simulation
# Possible values are: :sigma, :tau
imageplots = [:sigma]

# Whether to save GIFs of the simulation instead of individual image files
# Not recommended for long simulations (memory might fill up)
savegif = false

# How many MCSs must have passed before we try to save any output of the simulation
# Can be used to remove overhead of output functions
# Non-positive values disable this feature
outputperiod::Integer = 10

# How often in MCSs (+ 'outputperiod') to display information about the simulation on the console
# Non-positive values disable this feature
infoperiod::Integer = 1000

# How often in MCSs, (+ 'outputperiod') we save an image (or frame if 'savegif' is true) of the simulation
# Non-positive values disable this feature
imageperiod::Integer = 1000

# Framerate for the display updates (in real life time)
# Non-positive values disable this feature
# Low framerate may result in an unresponsive window
# This parameter is also affected by 'outputperiod'
displayframerate::Integer = 30

# How big is the window of the real-time display
displaysize::Integer = 350

# Whether to draw borders between cells
drawcellborders = true
