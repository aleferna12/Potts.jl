cellcolors = [RGB(rand(3)...) for _ in 1:100]

function plotsimulation!(cells::Cells, filename::AbstractString)
    img = fill(colorant"white", size(cells.matrix)...)
    # drawcells!(img, cells, colorant"gray")
    drawcells!(img, cells)
    drawcellborders!(img, cells, colorant"black")
    save(filename, img[end:-1:1, :])  # Reverse the image to comply with Images.jl standard
end

function drawcells!(img, cells, color::RGB)
    img[cells.matrix .> 0] .= color
end

function drawcells!(img, cells, colors=[cellcolors[attrset.sigma % length(cellcolors) + 1] for attrset in cells.attrsets])
    mask = cells.matrix .> 0
    img[mask] = colors[cells.matrix[mask]]
end

function drawcellborders!(img, cells, color::RGB)
    for edge in cells.edgeset
        for pos in [edge[1], edge[2]]
            sigma = cells.matrix[pos]
            if 1 ∉ [pos.x, pos.y] && sigma != 0 && (cells.matrix[pos.x - 1, pos.y] != sigma || cells.matrix[pos.x, pos.y - 1] != sigma)
                img[pos] = color
    end end end
end
