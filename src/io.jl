const SIGMACOLORS = [RGB(rand(3)...) for _ in 1:100]

function setupoutput(;infoperiod, imageperiod, imageplots, displayframerate, displaysize, kwargs...)
    outputobjs = (
        displayframerate > 0 ? makegui(length(imageplots), displaysize) : Dict(),
        IterationTimer(max(0, infoperiod), set=true, active=infoperiod > 0),
        IterationTimer(max(0, imageperiod), set=true, active=imageperiod > 0),
        Timer(Millisecond(max(0, round(1000 / displayframerate))), set=true, active=displayframerate > 0)
    )
    outputobjs
end

function output!(model::AbstractCPM,
                 time,
                 gui,
                 infocounter,
                 savecounter,
                 displaytimer)
    if fire!(infocounter, time)
        println("Timestep: $time")
    end

    savenow = fire!(savecounter, time)
    displaynow = fire!(displaytimer)
    if savenow || displaynow
        imgdict = simulationimages(getenv(model), model[:imageplots], model[:cellcolors], model[:drawcellborders], model[:drawcellcenters])
        if savenow
            dir = joinpath(model[:simdir], model[:imagesdirname])
            if model[:savegif]
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
            drawcellcenters!(img, livingcells(env), colorant"red")
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

drawsigmas!(img, env) = drawcells!(img, getmatrix(env), [SIGMACOLORS[sigma % length(SIGMACOLORS) + 1] for sigma in (getsigma(cell) for cell in livingcells(env))])
drawtaus!(img, env, taucolors::Vector) = drawcells!(img, getmatrix(env), [taucolors[tau] for tau in (gettau(cell) for cell in livingcells(env))])

function drawcellborders!(img, env::Environment, color::RGB)
    for edge in getedgeset(env)
        for pos in [edge[1], edge[2]]
            sigma = getmatrix(env)[pos]
            if 1 ∉ [pos.x, pos.y] && sigma != 0 && (getmatrix(env)[getx(pos) - 1, gety(pos)] != sigma || getmatrix(env)[getx(pos), gety(pos) - 1] != sigma)
                img[pos] = color
end end end end

function drawcellcenters!(img, cells, color)
    for center in (getcenter(cell) for cell in cells)
        img[round(Int, getx(center)), round(Int, gety(center))] = color
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

function endmessage(model::AbstractCPM, starttime)
    env = getenv(model)
    rnow = now()
    message = "
    Population:
        Total: $(length(livingcells(env)))
    "
    taus = (gettau(cell) for cell in livingcells(env))
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