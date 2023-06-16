abstract type AbstractCPM{ET <: Environment} end
Base.getindex(model::AbstractCPM, param::Symbol) = getparam(model, param)

function setup!(model::AbstractCPM)
    for tau in 1:model[:taus]
        cell_len = model[:cell_length][tau]
        for _ in 1:model[:ncells][tau]
            divtime = model[:divtime][tau]
            initcell!(getenv(model),
                      cell_len,
                      getrandompos(getmatrix(getenv(model)), ceil(Int, cell_len / 2)),
                      tau=tau, 
                      targetarea=model[:targetcellarea][tau],
                      divtimer=IterationTimer(divtime, active=divtime > 0))
end end end

function step!(model::AbstractCPM)
    updatehamiltonian!(model)
    updateattributes!(getenv(model), model[:targetcellarea], model[:divtargetcellarea])
    select!(getenv(model))
    reproduce!(getenv(model))
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
            copyspin!(env, edge) # TODO: this takes very long! Can it be optimized?
end end end

"Computes the local energy difference in case a lattice copy event (represented by an edge) occurs."
function deltahamiltonian(model::AbstractCPM, copyattempt::Edge)
    env = getenv(model)
    deltaH = 0.
    sigma1, sigma2 = getsigma(env, copyattempt[1]), getsigma(env, copyattempt[2])
    if sigma1 != MED
        cell = getcell(env, sigma1)
        deltaH += deltaHsize(getarea(cell), 1, gettargetarea(cell), model[:sizelambda])
    end
    if sigma2 != MED
        cell = getcell(env, sigma2)
        deltaH += deltaHsize(getarea(cell), -1, gettargetarea(cell), model[:sizelambda])
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
    if sigma1 == BORDER
        return sigma2 == MED ? 0. : model[:borderenergy]
    elseif sigma1 == MED
        return model[:adhesionmedium][gettau(getcell(getenv(model), sigma2))]
    end
    celladhenergy(getcell(getenv(model), sigma1), getcell(getenv(model), sigma2), model[:adhesiontable])
end

celladhenergy(cell1::AbstractCell, cell2::AbstractCell, adhesiontable) = adhesiontable[gettau(cell1), gettau(cell2)] * 2

function acceptcopy(deltaH, bolzmanntemp)
    if deltaH < 0
        return true
    end
    rand() < ℯ^(-deltaH / bolzmanntemp)
end

function copyspin!(env::Environment, edge::Edge)
    sigma1, sigma2 = getsigma(env, edge[1]), getsigma(env, edge[2])
    setsigma!(env, edge[2], sigma1)
    exchangepos!(env, sigma2, sigma1, edge)
    updateedges!(env, edge[2])
end

function exchangepos!(env::Environment, oldsigma, newsigma, edge)
    if newsigma > 0
        gainpos!(getcell(env, newsigma), edge[2])
    end
    if oldsigma > 0
        losepos!(getcell(env, oldsigma), edge[2])
    end
    nothing
end

"Updates any stale edges around 'pos'."
function updateedges!(env::Environment, pos::MatrixPos)
    sigmapos = getsigma(env, pos)
    for neigh in moore_neighbors(pos)
        sigmaneigh = getsigma(env, neigh)
        if sigmapos == sigmaneigh
            removeedges!(getedgeset(env), pos, neigh)
        elseif sigmaneigh != BORDER
            addedges!(getedgeset(env), pos, neigh)
end end end

struct CPM{ET, PT} <: AbstractCPM{ET}
    env::ET
    params::PT

    function CPM(env, params, typecheck::Bool=true)
        if typecheck
            warntype(typeof(env))
            warntype(typeof(params))
        end
        new{typeof(env), typeof(params)}(env, params)
    end
end
"Gets the environment associated with the model."
getenv(model::CPM) = model.env
"Gets all parameters stored in the model as a NamedTuple."
getparams(model::CPM) = model.params
getparam(model::CPM, param::Symbol) = getproperty(getparams(model), param)