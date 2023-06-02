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
function Genome(nreceptors::Integer, nreggenes::Integer, noutgenes::Integer, scale::Float64)
    recep = [Receptor(i, 0, scale) for i in 1:nreceptors]
    reg = [RegGene(i, -1, -1, 2 * rand() - 1) for i in 1:nreggenes]
    out = [OutGene(i, -1, -1, 2 * rand() - 1) for i in 1:noutgenes]
    Genome(recep, reg, out, [(2 * rand() - 1) for _ in range(1, nreggenes * (nreceptors + nreggenes + noutgenes))])
end
getreceptors(genome::Genome) = genome.receptors
getreggenes(genome::Genome) = genome.reggenes
getoutgenes(genome::Genome) = genome.outgenes

"Returns the weight of the influence of 'source' on the expression state of 'target'."
function getweight(genome::Genome, source::Receptor, target::RegGene)
    nrecep = length(getreceptors(genome))
    genome.weights[nrecep * (target.id - 1) + source.id]
end

function getweight(genome::Genome, source::RegGene, target::RegGene)
    nrecep, nreg = length(getreceptors(genome)), length(getreggenes(genome))
    genome.weights[nreg * nrecep + nreg * (target.id - 1) + source.id]
end

function getweight(genome::Genome, source::RegGene, target::OutGene)
    nrecep, nreg, nout = length(getreceptors(genome)), length(getreggenes(genome)), length(getoutgenes(genome))
    genome.weights[nreg * nrecep + nreg ^ 2 +  nout * (target.id - 1) + source.id]
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
    for gene in [getreggenes(genome); getoutgenes(genome)]
        gene.exprsig = gene.nextexprsig
    end
end