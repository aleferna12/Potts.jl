function setup()
    cells = Cells(PARAMS.fieldsize)
    for _ in 1:PARAMS.ncells
        spawncell!(cells, PARAMS.cell_length)
    end
    cells
end

function step(time::Integer, cells::Cells)
    updatehamiltonian!(cells)
    # Iterate the matrix and use the sigmas to calculate and update cell attributes
    # updateattributes!(cells)
end

function updatehamiltonian!(cells::Cells)
    # TODO: Check c++ code to understand how we used to do it
    # I think they are adding to a 'loop' variable but this variable is an int and they are adding floats (so doing nothing essentially)
    tovisit = length(cells.edgeset) / 8  # 8 because moore_neighbors gives 8 neighbours
    for _ in 1:floor(tovisit)
        edge = rand(cells.edgeset)
        if getsigma(cells, edge[1]) == getsigma(cells, edge[2])
            continue
        end

        dH = deltahamiltonian(cells, edge)
        if acceptcopy(dH, PARAMS.boltzmanntemp)
            copyspin!(cells, edge)
end end end

"Computes the local energy difference in case a lattice copy event (represented by an edge) occurs."
function deltahamiltonian(cells::Cells, copyattempt::Edge)
    deltaH = 0.
    sigma1, sigma2 = cells.matrix[[copyattempt...]]
    if sigma1 != 0
        deltaH += deltaHsize(getarea(cells, sigma1), 1, PARAMS.targetcellarea, PARAMS.sizelambda)
    end
    if sigma2 != 0
        deltaH += deltaHsize(getarea(cells, sigma2), -1, PARAMS.targetcellarea, PARAMS.sizelambda)
    end
    deltaH + interface_adhenergy(cells, copyattempt[2], sigma1, PARAMS.adhesiontable, PARAMS.adhesionmedium, PARAMS.borderenergy)
end

deltaHsize(area, deltaarea, targetarea, sizelambda) = sizelambda * deltaarea * (2 * (area - targetarea) + deltaarea)

function interface_adhenergy(cells::Cells, pos::MatrixPos, newsigma, adhesiontable, adhesionmedium, borderenergy)::Float64
    sigmapos = getsigma(cells, pos)
    energy = 0.
    for neigh in moore_neighbors(pos)
        sigmaneigh = getsigma(cells, neigh)
        if sigmaneigh < 0
            energy += (newsigma == 0 ? 0 : borderenergy) - (sigmapos == 0 ? 0 : borderenergy)
        else
            energy +=
            getadhesion(cells, newsigma, sigmaneigh, adhesiontable, adhesionmedium) -
            getadhesion(cells, sigmapos, sigmaneigh, adhesiontable, adhesionmedium)
    end end
    energy
end

function acceptcopy(deltaH, bolzmanntemp)
    if deltaH < 0
        return true
    end
    rand() < ℯ^(-deltaH / bolzmanntemp)
end

function copyspin!(cells::Cells, edge::Edge)
    sigma1, sigma2 = getsigma(cells, edge[1]), getsigma(cells, edge[2])
    setsigma(cells, edge[2], sigma1)
    if sigma1 != 0
        cells.attrsets[sigma1].area += 1
    end
    if sigma2 != 0
        cells.attrsets[sigma2].area -= 1
    end

    added_edges = Set{Edge}()
    for neigh in moore_neighbors(edge[2])
        sigmaneigh = getsigma(cells, neigh)
        if sigma1 == sigmaneigh
            removeedges!(cells.edgeset, edge[2], neigh)
        elseif sigmaneigh > -1
            addedges!(cells.edgeset, edge[2], neigh)
            addedges!(added_edges, edge[2], neigh)
    end end
    added_edges
end

function updateattributes!(cells::Cells)
    error("Not implemented")
end