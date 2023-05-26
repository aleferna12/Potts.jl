const SIGMACOLORS = [RGB(rand(3)...) for _ in 1:100]
const TAUCOLORS = [colorant"darkgray", colorant"lightgreen"]  # TODO: make parameter

function setupoutput(params)
    outputobjs = (
        params.displayframerate > 0 ? makegui(length(params.imageplots), params.displaysize) : Dict(),
        IterationTimer(max(0, params.infoperiod), set=true, active=params.infoperiod > 0),
        IterationTimer(max(0, params.imageperiod), set=true, active=params.imageperiod > 0),
        Timer(Millisecond(max(0, round(1000 / params.displayframerate))), set=true, active=params.displayframerate > 0)
    )
    outputobjs
end

function output!(env,
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
        imgdict = simulationimages(env, params.imageplots, params.cellcolors, params.drawcellborders, params.drawcellcenters)
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

function simulationimages(env::Environment, plots::Vector{Symbol}, cellcolors, drawcellborders, drawcellcenters)
    imgdict = Dict{Symbol, Matrix}()
    for plot in plots
        img = fill(colorant"white", size(getmatrix(env)))
        if plot === :sigma
            drawsigmas!(img, env)
        elseif plot === :tau
            drawtaus!(img, env, cellcolors)
        end
        if drawcellborders
            drawcellborders!(img, env, colorant"black")
        end
        if drawcellcenters
            drawcellcenters!(img, getcells(env), colorant"red")
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

function drawcells!(img, matrix, color::RGB)
    img[matrix .> 0] .= color
end
function drawcells!(img, matrix, colors::Vector)
    mask = matrix .> 0
    img[mask] = colors[matrix[mask]]
end

drawsigmas!(img, env) = drawcells!(img, getmatrix(env), [SIGMACOLORS[sigma % length(SIGMACOLORS) + 1] for sigma in getsigmas(getcells(env))])
drawtaus!(img, env, taucolors::Vector) = drawcells!(img, getmatrix(env), [taucolors[Int(attrset.tau)] for attrset in getcells(env).cellattrs])

function drawcellborders!(img, env::Environment, color::RGB)
    for edge in getedgeset(env)
        for pos in [edge[1], edge[2]]
            sigma = getmatrix(env)[pos]
            if 1 ∉ [pos.x, pos.y] && sigma != 0 && (getmatrix(env)[pos.x - 1, pos.y] != sigma || getmatrix(env)[pos.x, pos.y - 1] != sigma)
                img[pos] = color
end end end end

function drawcellcenters!(img, cells, color)
    for center in getcenters(cells)
        img[round(Int, center.x), round(Int, center.y)] = color
    end
end

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

function endmessage(env::Environment, starttime)
    rnow = now()
    message = "
    Population:
        Total: $(length(getcells(env)))
    "
    taus = gettaus(getcells(env))
    for tau in unique(taus)
        message *= "\tTau $tau: $(count(==(tau), taus))\n"
    end

    message *= "
    Total run time: $(Time(0) + (rnow - starttime))
    Start time: $starttime
    End time: $rnow
    "
    print(replace(message, "\n    " => "\n"))
end