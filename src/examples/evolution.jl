const PARAMS = merge(DEFAULTPARAMS, @typednamedtuple((mu=0.1, mustd=0.05)))

Base.@kwdef mutable struct EvolvingCell <: EvolvableCell
    sigma::Int
    tau::Int
    area::Int
    targetarea::Float64
    center::Pos{Float64}
    alive::Bool = true
    divtimer::IterationTimer
    food::Float64 = 0
    genome::Genome = Genome(0, 2, 1, [])
end
getgenome(cell::EvolvingCell) = cell.genome
getfitness(cell::EvolvingCell) = getexpressionsignal(getoutgenes(getgenome(cell))[1])

function exampleevolution()
    run(CPM(Dish(EvolvingCell, PARAMS.fieldsize), reassign(PARAMS, replaceprevsim=true, divtime=[5000], ncells=[20], imageplots=[:evoltarget, :sigma])))
end

function select!(env::Dish{EvolvingCell})
    killed = @invoke select!(env::Environment)
    nalive = length(collect(livingcells(env)))
    if nalive > 50
        sortedfit = collect(sortedfitness(env))
        for cell in sortedfit[1:(nalive - 50)]
            kill!(env, cell)
            push!(killed, cell)
    end end
    killed
end

function simulationimages(env::Dish{EvolvingCell}, plots::Vector{Symbol}, cellcolors, drawcellborders, drawcellcenters)
    imgdict = @invoke simulationimages(env::Environment,
                                       filter(x -> x !== :evoltarget, plots),
                                       cellcolors, 
                                       drawcellborders, 
                                       drawcellcenters)
    if :evoltarget in plots
        img = fill(colorant"white", size(getmatrix(env)))
        drawevoltarget!(img, env)
        imgdict[:evoltarget] = img[end:-1:1, :]
    end
    imgdict
end

function drawevoltarget!(img, env)
    colors = [RGB(min(1, max(0, getfitness(cell) / 3))) for cell in getcells(env)]
    drawcells!(img, env, colors)
end

function maxfitness(env)
    maximum((getfitness(cell) for cell in livingcells(env)))
end

function sortedfitness(env)
    cellfit = [(cell, getfitness(cell)) for cell in livingcells(env)]
    sort!(cellfit, by=x -> x[2])
    (tup[1] for tup in cellfit)
end