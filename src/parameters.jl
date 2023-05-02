Base.@kwdef struct Parameters
    "Size of the field."
    fieldsize::Int = 100
    "Number of starting cells."
    ncells::Int = 15
    "Starting cell length."
    cell_length::Int = 5
    "Target area for cells."
    targetcellarea::Int = 25
    "Size constraint scale (how strongly cells deviating from 'targetcellarea' are penalized)."
    sizelambda::Float64 = 4
    "Symmetric table used to determine cell-cell adhesion energy. The indices are the taus of cells."
    adhesiontable::Matrix{Float64} = fill(16, 2, 2)
    "Adhesion energy of a cell with the medium."
    adhesionmedium::Float64 = 8
    "Energy on the interface with a border site."
    borderenergy::Float64 = 100
    "Temperature controling the 'thermodynamic' behaviour of the system."
    boltzmanntemp::Float64 = 12
    "How many MCSs to run the simulation for."
    endsim::Int = 1000
    "How often in MCSs we save an image of the sim."
    imageperiod::Int = 200
end

# This should be a const but Revise doesnt like consts so im keeping it like that for now 
PARAMS::Parameters = Parameters(  # TODO: Replace with readparameters function
    boltzmanntemp=12,
    targetcellarea=150,
    sizelambda=100000
)