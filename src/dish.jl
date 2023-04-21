function setup(;params...)
    cells = Cells(params[:size])
    for _ in 1:params[:ncells]
        spawncell!(cells, params[:cell_length])
    end
    cells
end

function step(time::Integer, cells::Cells; params...)
    updatehamiltonian!(cells; params...)
    # Iterate the matrix and use the sigmas to calculate and update cell attributes
    # updateattributes!(cells)
end

function updatehamiltonian!(cells::Cells; params...)
    edgeset = copy(cells.edgeset)
    for edge in (pop!(edgeset) for _ in 1:length(edgeset))
        dH = deltahamiltonian(edge, cells; params...)

        if !acceptcopy(dH, params[:boltzmanntemp])
            continue
        end

        sigmas = cells.matrix[edge[1]], cells.matrix[edge[2]]
        cells.matrix[edge[2]] = sigmas[1]
        if sigmas[1] != 0
            cells.attrsets[sigmas[1]].area -= 1
        end
        if sigmas[2] != 0
            cells.attrsets[sigmas[2]].area += 1
        end
        removeedges!(cells.edgeset, edge...)
        for neigh in moore_neighbors(edge[2])
            if cells.matrix[neigh] != cells.matrix[edge[2]]
                addedges!(cells.edgeset, edge[2], neigh)
end end end end

"Computes the local energy difference in case a lattice copy event (represented by an edge) occurs."
function deltahamiltonian(copyattempt::Edge, cells; params...)
    sigmas = orderpair(cells.matrix[copyattempt[1]], cells.matrix[copyattempt[2]])
    # TODO Fix this should take into account order size2 - size1, not bigger only
    if sigmas[1] == 0
        dHsize = deltasize(cells.attrsets[sigmas[2]].area, params[:targetarea], params[:sizelambda])
        return float(dHsize + params[:adhesionmedium])
    end
    
    attrs1, attrs2 = cells.attrsets[sigmas[1]], cells.attrsets[sigmas[2]]
    dHsize = deltasize(attrs1.area, attrs2.area, params[:targetarea], params[:sizelambda])
    dHadh = deltaadhesion(attrs1.tau, attrs2.tau, params[:adhesiontable])
    float(dHsize + dHadh)
end

deltasize(area, targetarea, lambda::Real) = lambda * (area - targetarea)^2
deltasize(area1, area2, targetarea, lambda::Real) = lambda * ((area2 - targetarea)^2 - (area1 - targetarea)^2)
deltaadhesion(tau1, tau2, adhesiontable::Matrix) = adhesiontable[Integer(tau1), Integer(tau2)]

function acceptcopy(deltaH, bolzmanntemp)
    if deltaH < 0
        return true
    end
    rand() < ℯ^(-deltaH / bolzmanntemp)
end

function updateattributes!(cells::Cells)
    error("Not implemented")
end