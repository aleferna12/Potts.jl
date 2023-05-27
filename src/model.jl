abstract type AbstractCPM end
Base.getindex(model::AbstractCPM, param::Symbol) = getparam(model, param)

struct CPM{ET, PT} <: AbstractCPM where ET <: Environment
    env::ET
    params::PT

    function CPM(env, params, typecheck::Bool=true)
        if typecheck
            if !fullyconcrete(typeof(env))
                @warn "'$env' may contain abstract fields in it's type hierarchy"
            end
            if !fullyconcrete(typeof(params))
                @warn "'$params' may contain abstract fields in it's type hierarchy"
            end
        end
        new{typeof(env), typeof(params)}(env, params)
    end
end
"Gets the environment associated with the model."
getenv(model::CPM) = model.env
"Gets all parameters stored in the model as a NamedTuple."
getparams(model::CPM) = model.params
getparam(model::CPM, param::Symbol) = getproperty(getparams(model), param)

function setup!(model::AbstractCPM)
    for tau in 1:model[:taus]
        for _ in 1:model[:ncells][tau]
            spawncell!(getenv(model), tau, model[:cell_length][tau])
end end end

function step!(model::AbstractCPM, time::Integer)
    updatehamiltonian!(model)
    # Iterate the matrix and use the sigmas to calculate and update cell attributes
    # updateattributes!(cells)
end

function updatehamiltonian!(model::AbstractCPM)
    # TODO: Check c++ code to understand how we used to do it
    # I think they are adding to a 'loop' variable but this variable is an int and they are adding floats (so doing nothing essentially)
    env = getenv(model)
    tovisit = length(getedgeset(env)) / 8  # 8 because moore_neighbors gives 8 neighbours
    for _ in 1:ceil(tovisit)
        edge = rand(getedgeset(env)) # TODO: this takes very long! Can it be optimized?
        if getsigma(env, edge[1]) == getsigma(env, edge[2])
            continue
        end

        dH = deltahamiltonian(model, edge)
        if acceptcopy(dH, model[:boltzmanntemp])
            copyspin!(env, edge)
end end end

"Computes the local energy difference in case a lattice copy event (represented by an edge) occurs."
function deltahamiltonian(model::AbstractCPM, copyattempt::Edge)
    env = getenv(model)
    deltaH = 0.
    sigma1, sigma2 = getsigma(env, copyattempt[1]), getsigma(env, copyattempt[2])
    if sigma1 != 0
        deltaH += deltaHsize(getarea(getcell(env, sigma1)), 1, model[:targetcellarea][gettau(getcell(env, sigma1))], model[:sizelambda])
    end
    if sigma2 != 0
        deltaH += deltaHsize(getarea(getcell(env, sigma2)), -1, model[:targetcellarea][gettau(getcell(env, sigma2))], model[:sizelambda])
    end
    neighsigmas = [getsigma(env, nsig) for nsig in moore_neighbors(copyattempt[2])]
    deltaH + deltaHadhenergy(model, sigma2, sigma1, neighsigmas)
end

deltaHsize(area, deltaarea, targetarea, sizelambda) = sizelambda * deltaarea * (2 * (area - targetarea) + deltaarea)

function deltaHadhenergy(model::AbstractCPM,
                         oldsigma,
                         newsigma,
                         neighsigmas)
    energy = 0.
    for neighsigma in neighsigmas
        energy += getadhenergy(model, newsigma, neighsigma)
        energy -= getadhenergy(model, oldsigma, neighsigma)
    end
    energy
end

function getadhenergy(model::AbstractCPM, sigma1::Integer, sigma2::Integer)::Float64
    if sigma1 == sigma2
        return 0.
    end
    sigma1, sigma2 = ordered(sigma1, sigma2)
    if sigma1 < 0
        return sigma2 == 0 ? 0. : model[:borderenergy]
    elseif sigma1 == 0
        return model[:adhesionmedium][gettau(getcell(getenv(model), sigma2))]
    end
    model[:adhesiontable][gettau(getcell(getenv(model), sigma1)), gettau(getcell(getenv(model), sigma2))] * 2
end

function acceptcopy(deltaH, bolzmanntemp)
    if deltaH < 0
        return true
    end
    rand() < ℯ^(-deltaH / bolzmanntemp)
end

function copyspin!(env::Environment, edge::Edge)
    sigma1, sigma2 = getsigma(env, edge[1]), getsigma(env, edge[2])
    setsigma(env, edge[2], sigma1)
    updateattributes!(env, sigma2, sigma1, edge)
    updateedges!(env, sigma1, edge[2])
end