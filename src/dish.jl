function setup(params)
    cells = Cells(params.fieldsize)
    for _ in 1:params.ncells
        spawncell!(cells, params.cell_length)
    end
    cells
end

function step(cells::Cells, time::Integer, params)
    updatehamiltonian!(cells, params)
    # Iterate the matrix and use the sigmas to calculate and update cell attributes
    # updateattributes!(cells)
end

function updatehamiltonian!(cells::Cells, params)
    # TODO: Check c++ code to understand how we used to do it
    # I think they are adding to a 'loop' variable but this variable is an int and they are adding floats (so doing nothing essentially)
    tovisit = length(cells.edgeset) / 8  # 8 because moore_neighbors gives 8 neighbours
    for _ in 1:ceil(tovisit)
        edge = rand(cells.edgeset)
        if getsigma(cells, edge[1]) == getsigma(cells, edge[2])
            continue
        end

        dH = deltahamiltonian(cells, edge, params)
        if acceptcopy(dH, params.boltzmanntemp)
            copyspin!(cells, edge)
end end end

"Computes the local energy difference in case a lattice copy event (represented by an edge) occurs."
function deltahamiltonian(cells::Cells, copyattempt::Edge, params)
    deltaH = 0.
    sigma1, sigma2 = cells.matrix[[copyattempt...]]
    if sigma1 != 0
        deltaH += deltaHsize(getarea(cells, sigma1), 1, params.targetcellarea, params.sizelambda)
    end
    if sigma2 != 0
        deltaH += deltaHsize(getarea(cells, sigma2), -1, params.targetcellarea, params.sizelambda)
    end
    deltaH + deltaHadhenergy(sigma2, sigma1, copyattempt[2], cells, params.adhesiontable, params.adhesionmedium, params.borderenergy)
end

deltaHsize(area, deltaarea, targetarea, sizelambda) = sizelambda * deltaarea * (2 * (area - targetarea) + deltaarea)

function deltaHadhenergy(oldsigma,
                         newsigma,
                         pos::MatrixPos,
                         cells,
                         adhargs...)
    energy = 0.
    for neighsigma in (getsigma(cells, npos) for npos in moore_neighbors(pos))
        energy += getadhenergy(cells, newsigma, neighsigma, adhargs...)
        energy -= getadhenergy(cells, oldsigma, neighsigma, adhargs...)
    end
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