const SIGMACOLORS = [RGB(rand(3)...) for _ in 1:100]
const TAUCOLORS = [colorant"darkgray", colorant"lightgreen"]  # TODO: make parameter

function setupoutput(params)
    outputobjs = (
        params.displayframerate > 0 ? makegui(length(params.imageplots), params.displaysize) : Dict(),
        Counter(max(0, params.infoperiod), cock=true, active=params.infoperiod > 0),
        Counter(max(0, params.imageperiod), cock=true, active=params.imageperiod > 0),
        Timer(Millisecond(max(0, round(1000 / params.displayframerate))), cock=true, active=params.displayframerate > 0)
    )
    outputobjs
end

function output!(cells,
                 time,
                 gui,
                 infocounter,
                 savecounter,
                 displaytimer,
                 params)
    if fire!(infocounter, time)
        println("Timestep: $time")
    end

    savenow = fire!(savecounter, time)
    displaynow = fire!(displaytimer)
    if savenow || displaynow
        imgdict = simulationimages(cells, params.imageplots, params.drawcellborders)
        if savenow
            dir = joinpath(params.simdir, params.imagesdirname)
            if params.savegif
                save_simulationimages(imgdict, dir)
            else
                save_simulationimages(imgdict, dir, time)
        end end
        if displaynow
            display_simulationimages(imgdict, gui)
    end end
end

function simulationimages(cells::Cells, plots::Vector{Symbol}, drawcellborders)
    imgdict = Dict{Symbol, Matrix}()
    for plot in plots
        img = fill(colorant"white", size(cells.matrix))
        if plot == :sigma
            drawsigmas!(img, cells)
        elseif plot == :tau
            drawtaus!(img, cells, TAUCOLORS)
        end
        if drawcellborders
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
    sleep(0.000001) # This is the time frame in which the user can interact with the display
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

function makegui(nplots::Integer, canvassize)
    rows = nplots == 1 ? 1 : 2
    cols = ceil(Int, nplots / 2)
    gui = imshow_gui((canvassize, canvassize), (rows, cols))
    Gtk.resize!(gui["window"], canvassize * cols, canvassize * rows)
    Gtk.showall(gui["window"])
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