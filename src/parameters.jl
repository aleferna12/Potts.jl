Base.@kwdef struct Parameters
    ncells::Int = 50
    fieldsize::Int = 100
    cell_length::Int = 5
    targetcellarea::Int = 25
    sizelambda::Float64 = 4
    adhesiontable::Matrix{Float64} = zeros(2, 2)
    adhesionmedium::Float64 = 0
    boltzmanntemp::Float64 = 16
end

const PARAMS = Parameters()