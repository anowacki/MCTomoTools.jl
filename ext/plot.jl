function Makie.plot(
    chains::Union{Chain,<:AbstractArray{<:Chain}},
    field::Symbol;
    figure=(),
    axis=(),
    kwargs...
)
    fig = Makie.Figure(; figure...)
    ax, pl = Makie.plot(fig[1,1], chains, field; axis, kwargs...)
    Makie.FigureAxisPlot(fig, ax, pl)
end

function Makie.plot(
    gp::Makie.GridPosition,
    chains::Union{Chain,<:AbstractArray{<:Chain}},
    field::Symbol;
    axis=(),
    kwargs...
)
    axis = Makie.Axis(gp;
        xlabel="Step number",
        ylabel=uppercasefirst(String(field)),
        axis...
    )
    pl = Makie.plot!(axis, chains, field; kwargs...)
    Makie.AxisPlot(axis, pl)
end

function Makie.plot!(
    axis::Makie.Axis,
    chains::Union{Chain,<:AbstractArray{<:Chain}},
    field::Symbol;
    thin=1000,
    burnin=0,
    colormap=:turbo,
    series=(),
)
    chains isa Chain && (chains = [chains])
    nchains = length(chains)

    isempty(chains) &&
        throw(ArgumentError("chains is empty; must pass at least one chain in"))

    # Extract all samples
    isample1 = burnin + 1
    isample2 = maximum(length, chains)

    # Ignore `thin` if it's too large
    nsamples = isample2 - isample1 + 1
    thin = if thin > nsamples/100
        ceil(Int, nsamples/100)
    else
        thin
    end

    # Thin out the data and chop to the desired window
    samples = [@view chain.samples[max(burnin + 1, isample1):thin:min(end, isample2)]
        for chain in chains]
    indices = [eachindex(chain.samples)[max(burnin + 1, isample1):thin:min(end, isample2)]
        for chain in chains]

    # Trace data
    data = if field in fieldnames(MCTomoTools.RawSample)
        [getproperty.(s, field) for s in samples]
    elseif field === :loglikelihood
        # Need to ensure everything is positive
        likelihoods = [-getproperty.(s, :likelihood) for s in samples]
        for ll in likelihoods
            for (i, val) in pairs(ll)
                ll[i] = val < 0 ? NaN : log(val)
            end
        end
        likelihoods
    else
        throw(ArgumentError("unknown field name ':$field'"))
    end

    # Trace annotations, on right of last points
    annot_xs = last.(indices)
    annot_ys = last.(data)
    annot_text = " ".*string.(1:nchains)

    # Plot traces as a series so this is faster for multiple traces
    points = [Makie.Point2f.(x, y) for (x, y) in zip(indices, data)]
    # FIXME: Ugly internal field access to get colours
    colors = Makie.cgrad(
        colormap,
        range(0, 1, length=(nchains + 1));
        categorical=true
    ).colors.colors

    # Trace plot
    pl = Makie.series!(axis, points; labels=string.(1:nchains), color=colors, series...)

    # Annotations
    Makie.text!(axis, annot_xs, annot_ys;
        text=annot_text, color=colors, align=(:left, :baseline)
    )

    pl
end
