const SIGMACOLORS = [RGB(rand(3)...) for _ in 1:100]
const TAUCOLORS = [colorant"darkgray", colorant"lightgreen"]  # TODO: make parameter

function output(cells, gui, time)
    if PARAMS.infoperiod > 0 && time % PARAMS.infoperiod == 0
        println("Timestep: $time")
    end

    savetime = PARAMS.imageperiod > 0 && time % PARAMS.imageperiod == 0
    displaytime = PARAMS.displayperiod > 0 && time % PARAMS.displayperiod == 0
    if savetime || displaytime
        imgdict = simulationimages(cells, PARAMS.imageplots)
        if savetime
            dir = joinpath(PARAMS.simdir, PARAMS.imagesdirname)
            if PARAMS.savegif
                save_simulationimages(imgdict, dir)
            else
                save_simulationimages(imgdict, dir, time)
        end end
        if displaytime
            display_simulationimages(imgdict, gui)
end end end

function simulationimages(cells::Cells, plots::Vector{Symbol})
    imgdict = Dict{Symbol, Matrix}()
    for plot in plots
        img = fill(colorant"white", size(cells.matrix))
        if plot == :sigma
            drawsigmas!(img, cells)
        elseif plot == :tau
            drawtaus!(img, cells, TAUCOLORS)
        end
        if PARAMS.drawcellborders
            drawcellborders!(img, cells, colorant"black")
        end
        imgdict[plot] = img[end:-1:1, :]
    end
    imgdict
end

function save_simulationimages(imgdict, dir)
    for (plot, img) in imgdict
        filepath = joinpath(dir, "$plot.gif")
        if !isfile(filepath)
            save(filepath, img)
        else
            image3d = load(filepath)
            save(filepath, cat(image3d, img, dims=3))
end end end
function save_simulationimages(imgdict, dir, time)
    for (plot, img) in imgdict
        filename = joinpath(dir, "$plot$time.png")
        save(filename, img)
end end

function display_simulationimages(imgdict, gui)
    if length(imgdict) == 1
        imshow(gui["canvas"], values(imgdict) |> first)
    else
        for (i, img) in values(imgdict) |> enumerate
            imshow(gui["canvas"][i], img)
    end end
    sleep(0.0001) # Forces the UI to update 
end

function drawcells!(img, cells, color::RGB)
    img[cells.matrix .> 0] .= color
end
function drawcells!(img, cells, colors::Vector)
    mask = cells.matrix .> 0
    img[mask] = colors[cells.matrix[mask]]
end

drawsigmas!(img, cells) = drawcells!(img, cells, [SIGMACOLORS[attrset.sigma % length(SIGMACOLORS) + 1] for attrset in cells.attrsets])
drawtaus!(img, cells, taucolors::Vector) = drawcells!(img, cells, [taucolors[Int(attrset.tau)] for attrset in cells.attrsets])

function drawcellborders!(img, cells, color::RGB)
    for edge in cells.edgeset
        for pos in [edge[1], edge[2]]
            sigma = cells.matrix[pos]
            if 1 ∉ [pos.x, pos.y] && sigma != 0 && (cells.matrix[pos.x - 1, pos.y] != sigma || cells.matrix[pos.x, pos.y - 1] != sigma)
                img[pos] = color
end end end end

function makegui(displayperiod::Integer, nplots::Integer, canvassize)
    gui::Dict{String, Any} = Dict()
    if displayperiod > 0
        rows = nplots == 1 ? 1 : 2
        cols = ceil(Int, nplots / 2)
        gui = imshow_gui((canvassize, canvassize), (rows, cols))
        Gtk.resize!(gui["window"], canvassize * cols, canvassize * rows)
        Gtk.showall(gui["window"])
    end
    gui
end

function makesimdirs(simdir, subdirs, replaceprevsim)
    if isdir(simdir)
        if replaceprevsim
            rm(simdir, recursive=true)
        else
            error("directory $simdir already exists")
    end end
    for subdir in subdirs
        joinpath(simdir, subdir) |> mkpath
end end