function setup()
    cells = Cells(PARAMS.fieldsize)
    for _ in 1:PARAMS.ncells
        spawncell!(cells, PARAMS.cell_length)
    end
    calculate_adhesion_matrix!(cells, PARAMS.adhesiontable, PARAMS.adhesionmedium)
    cells
end

function calculate_adhesion_matrix!(cells::Cells, adhesiontable::Matrix, adhesionmedium::AbstractFloat)
    for edge in cells.edgeset
        cells.adhmatrix[edge[1]] = interface_adhenergy(cells, edge[1], adhesiontable, adhesionmedium)
        cells.adhmatrix[edge[2]] = interface_adhenergy(cells, edge[2], adhesiontable, adhesionmedium)
end end

function step(time::Integer, cells::Cells)
    updatehamiltonian!(cells)
    # Iterate the matrix and use the sigmas to calculate and update cell attributes
    # updateattributes!(cells)
end

function updatehamiltonian!(cells::Cells)
    edgeset = copy(cells.edgeset)
    for edge in (pop!(edgeset) for _ in 1:length(edgeset))
        dH, new_adhenergy = deltahamiltonian(cells, edge)

        if acceptcopy(dH, PARAMS.boltzmanntemp)
            copyspin!(cells, edge, new_adhenergy)
end end end

"Computes the local energy difference in case a lattice copy event (represented by an edge) occurs."
function deltahamiltonian(cells::Cells, copyattempt::Edge)
    if copyattempt[1].x < copyattempt[2].x
        return -1.
    else
        return Inf
    end
    deltaH = 0.
    sigma1, sigma2 = cells.matrix[[copyattempt...]]
    if sigma1 != 0
        deltaH += deltaHsize(getarea(cells, sigma1), -1, PARAMS.targetcellarea, PARAMS.sizelambda)
    end
    if sigma2 != 0
        deltaH += deltaHsize(getarea(cells, sigma2), 1, PARAMS.targetcellarea, PARAMS.sizelambda)
    end
    new_adhenergy = interface_adhenergy(cells, copyattempt[2], PARAMS.adhesiontable, PARAMS.adhesionmedium, sigma2)
    deltaH += new_adhenergy - getadhenergy(cells, copyattempt[2])
    deltaH, new_adhenergy
end

deltaHsize(area, deltaarea, targetarea, sizelambda) = sizelambda * deltaarea * (2 * area - 2 * targetarea + deltaarea)

function interface_adhenergy(cells::Cells, pos::MatrixPos, adhesiontable, adhesionmedium, sigmapos=getsigma(cells, pos))::Float64
    neighs = getneighbors(cells, pos)
    if sigmapos == 0
        return adhesionmedium * length(neighs)
    end

    energy = 0.
    for neigh in neighs
        energy += getadhesion(cells, sigmapos, getsigma(cells, neigh), adhesiontable, adhesionmedium)
    end
    energy
end

function acceptcopy(deltaH, bolzmanntemp)
    if deltaH < 0
        return true
    end
    rand() < ℯ^(-deltaH / bolzmanntemp)
end

function copyspin!(cells::Cells, edge::Edge, new_adhenergy)
    sigma1, sigma2 = getsigma(cells, edge[1]), getsigma(cells, edge[2])
    setsigma(cells, edge[2], sigma1)
    setadhenergy(cells, edge[2], new_adhenergy)
    if sigma1 != 0
        cells.attrsets[sigma1].area -= 1
    end
    if sigma2 != 0
        cells.attrsets[sigma2].area += 1
    end

    removeedges!(cells.edgeset, edge...)
    for neigh in getneighbors(cells, edge[2])
        if sigma2 != getsigma(cells, neigh)
            addedges!(cells.edgeset, edge[2], neigh)
end end end

function updateattributes!(cells::Cells)
    error("Not implemented")
end