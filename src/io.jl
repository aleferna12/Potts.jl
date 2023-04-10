using PlotlyJS

function plotsimulation!(cellfield::Matrix, filename::AbstractString)
    trace = heatmap(z=cellfield)
    layout = Layout(
        width=size(cellfield)[1],
        height=size(cellfield)[2],
        yaxis=attr(scaleanchor="x", scaleratio=1, constrain="domain"),
        xaxis=attr(constrain="domain")
    )
    p = plot(trace, layout)

    savefig(p, filename)
end