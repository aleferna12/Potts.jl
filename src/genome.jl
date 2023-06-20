mutable struct Receptor
    id::Int
    value::Float64
    scale::Float64
end
Base.read(receptor::Receptor) = receptor.value * receptor.scale
write!(receptor::Receptor, val) = receptor.value = val

abstract type AbstractGene end

mutable struct RegGene
    id::Int
    exprsig::Float64
    nextexprsig::Float64
    threshold::Float64
end
getexpressionsignal(gene::RegGene) = gene.exprsig
getthreshold(gene::RegGene) = gene.threshold
isexpressed(gene::RegGene) = getexpressionsignal(gene) > getthreshold(gene)

mutable struct OutGene
    id::Int
    exprsig::Float64
    nextexprsig::Float64
    threshold::Float64
end
getexpressionsignal(gene::OutGene) = gene.exprsig
getthreshold(gene::OutGene) = gene.threshold
isexpressed(gene::OutGene) = getexpressionsignal(gene) > getthreshold(gene)

struct Genome
    receptors::Vector{Receptor}
    reggenes::Vector{RegGene}
    outgenes::Vector{OutGene}
    weights::Vector{Float64}
end
function Genome(nreceptors::Integer, nreggenes::Integer, noutgenes::Integer, recepscales)
    recep = [Receptor(i, 0, recepscales[i]) for i in 1:nreceptors]
    reg = [RegGene(i, 0, 0, 2 * rand() - 1) for i in 1:nreggenes]
    out = [OutGene(i, 0, 0, 2 * rand() - 1) for i in 1:noutgenes]
    Genome(recep, reg, out, [(2 * rand() - 1) for _ in range(1, nreggenes * (nreceptors + nreggenes + noutgenes))])
end
getreceptors(genome::Genome) = genome.receptors
getreggenes(genome::Genome) = genome.reggenes
getoutgenes(genome::Genome) = genome.outgenes

function indexweight(genome::Genome, source::Receptor, target::RegGene)
    length(getreceptors(genome)) * (target.id - 1) + source.id
end

function indexweight(genome::Genome, source::RegGene, target::RegGene)
    nreg = length(getreggenes(genome))
    nreg * length(getreceptors(genome)) + nreg * (target.id - 1) + source.id
end

function indexweight(genome::Genome, source::RegGene, target::OutGene)
    nreg = length(getreggenes(genome))
    nreg * length(getreceptors(genome)) + nreg ^ 2 +  length(getoutgenes(genome)) * (target.id - 1) + source.id
end

"Returns the weight of the influence of 'source' on the expression state of 'target'."
function getweight(genome::Genome, source, target)
    genome.weights[indexweight(genome, source, target)]
end

function setweight!(genome::Genome, source, target, val)
    genome.weights[indexweight(genome, source, target)] = val
end

function reset!(genome::Genome)
    for reggene in getreggenes(genome)
        reggene.exprsig = 0
        reggene.nextexprsig = 0
    end
    for outgene in getoutgenes(genome)
        outgene.exprsig = 0
        outgene.nextexprsig = 0
    end
end

function update!(genome::Genome)
    updatereggenes!(genome)
    updateoutgenes!(genome)
    finishupdate!(genome)
end

function updatereggenes!(genome::Genome)
    for reg in getreggenes(genome)
        exprsig = 0.
        for recep in getreceptors(genome)
            exprsig += read(recep) * getweight(genome, recep, reg)
        end
        for reg2 in getreggenes(genome)
            isexpressed(reg2) && (exprsig += getweight(genome, reg2, reg))
        end
        reg.nextexprsig = exprsig
end end

function updateoutgenes!(genome::Genome)
    for out in getoutgenes(genome)
        exprsig = 0.
        for reg in getreggenes(genome)
            isexpressed(reg) && (exprsig += getweight(genome, reg, out))
        end
        out.nextexprsig = exprsig
end end

function finishupdate!(genome::Genome)
    for gene in getreggenes(genome)
        gene.exprsig = gene.nextexprsig
    end
    for gene in getoutgenes(genome)
        gene.exprsig = gene.nextexprsig
    end
end

function mutate!(genome::Genome, mu, mustd)
    dist = Normal(0, mustd)
    for reggene in getreggenes(genome)
        for recep in getreceptors(genome)
            if rand() < mu
                setweight!(genome, recep, reggene, getweight(genome, recep, reggene) + rand(dist))
            end
        end
        for reggene2 in getreggenes(genome)
            if rand() < mu
                setweight!(genome, reggene, reggene2, getweight(genome, reggene, reggene2) + rand(dist))
            end
        end
        if rand() < mu
            reggene.threshold += rand(dist)
        end
    end
    for outgene in getoutgenes(genome)
        for reggene in getreggenes(genome)
            if rand() < mu
                setweight!(genome, reggene, outgene, getweight(genome, reggene, outgene) + rand(dist))
            end
        end
        if rand() < mu
            outgene.threshold += rand(dist)
        end
    end
end