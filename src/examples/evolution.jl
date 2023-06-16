Base.@kwdef mutable struct EvolvingCell <: EvolvableCell
    sigma::Int
    tau::Int
    area::Int
    targetarea::Float64
    center::Pos{Float64}
    alive::Bool = true
    divtimer::IterationTimer
    food::Float64 = 0
    genome::Genome = Genome(0, 3, 1, ())
end

function examplehungry()
    params = defaultparams(CPM)
    run(CPM(Dish(EvolvingCell, params.fieldsize), reassign(params, replaceprevsim=true)))
end