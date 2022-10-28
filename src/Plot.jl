"""
# `Plot`

Module defining plotting recipes and helper functions for plotting.
"""
module Plot

using RecipesBase: RecipesBase, @recipe, @series, @userplot
using MCTomoTools: MCTomoTools, Chain, Grid, sample_average_grid

"""
    plot(chains::Vector{Chain}, field::Symbol; burnin=0, thin=100) -> ::Plot

Plot the value of `field` at every `thin`th sample in each chain in `chains`,
discarding the first `burnin` samples.

`chains` may be either a set of `Chain`s, or a single one.
"""
function plot end

@recipe function f(chains::Union{Chain,<:AbstractVector{<:Chain}}, field::Symbol;
    thin=1000, burnin=0
)
    chains isa Chain && (chains = [chains])
    nchains = length(chains)
    isempty(chains) &&
        throw(ArgumentError("chains is empty; must pass at least one chain in"))

    # If we have requested only a certain part of the plot, we only need
    # to extract that information
    isample1, isample2 = get(plotattributes, :xlims,
        (burnin + 1, maximum(length, chains)))

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

    # Trace labels
    label_xs = last.(indices)
    label_ys = last.(data)
    labels = [(string(" ", chain), 8, :left) for chain in chains]
    @show label_xs, label_ys, labels

    # Nicer defaults
    fontfamily --> "Helvetica"
    framestyle --> :box
    grid --> false
    legend --> false

    @series begin
        xguide --> "Sample number"
        yguide --> uppercasefirst(String(field))
        annotations := #=nchains > 1=# false ?
            [(label_x, label_y, label)
                for (label_x, label_y, label) in zip(label_xs, label_ys, labels)] :
            []
        indices, data
    end
end

"""
    plot(chain::Chain, grid::Grid, property::Symbol; x, y, z, burnin=0, thin=1000, method=:brute) -> ::Plot

Plot a mean cross-section through `grid` at either x=`x`, y=`y` or z=`z`,
for `property`, which may be `:vp`, `:vs` or `:density`.

`method` is passed to [`sample_average_grid`](@ref) and controls how
the cross-section is created.
"""
plot_section

@userplot Plot_Section

@recipe function f(gridsection::Plot_Section;
    # chain::Chain, grid::Grid, property::Symbol;
    burnin=0, thin=1000, x=nothing, y=nothing, z=nothing, method=:brute
)
    nargs = length(gridsection.args)
    3 <= nargs <= 5 ||
        throw(ArgumentError("arguments to plot_section are chain, grid, property[, sources[, receoivers]]"))
    chains, grid, property = gridsection.args
    (chains isa Chain || chains isa AbstractArray{<:Chain}) ||
        throw(ArgumentError("chain must of type Chain or Vector{Chain}"))
    grid isa Grid || throw(ArgumentError("chain must of type Grid"))
    property isa Symbol || throw(ArgumentError("chain must of type Symbol"))

    count(isnothing, (x, y, z)) == 2 ||
        throw(ArgumentError("one and only one of x, y or z must be given"))

    # Optional positional arguments
    sources = nargs >= 4 ? gridsection.args[4] : nothing
    receivers = nargs == 5 ? gridsection.args[5] : nothing

    # Actually make a new grid from the grid given in
    slice_dim = findfirst(!isnothing, (x, y, z))
    slice_name = (:x, :y, :z)[slice_dim]
    slice = if slice_name === :x
        Grid(x:1:x, grid.y, grid.z)
    elseif slice_name === :y
        Grid(grid.x, y:1:y, grid.z)
    else
        Grid(grid.x, grid.y, z:1:z)
    end

    # Now calculate mean and stdev grid
    μ, σ = sample_average_grid(chains, slice, property; burnin, thin, method)
    μ = dropdims(μ, dims=slice_dim)
    σ = dropdims(σ, dims=slice_dim)

    # Work out correct axes and things for this slice
    xcoords = slice_name === :x ? grid.y : grid.x
    ycoords = slice_name === :z ? grid.y : grid.z
    xlabel = (slice_name === :x ? "Northing / km" : "Easting / km")
    ylabel = (slice_name === :z ? "Northing / km" : "Depth / km")

    property_name = uppercasefirst(String(property))
    units = property === :density ? "g/cm^3" : "km/s"

    # Sensible plotting defaults
    layout --> (1, 2)
    fontfamily --> "Helvetica"
    framestyle --> :box
    grid --> false

    xlims --> extrema(xcoords)
    ylims --> extrema(ycoords)
    yflip --> slice_name !== :z

    # Cross sections: fix labels, aspect ratio, etc., at the end
    # Mean section
    @series begin
        subplot := 1
        seriestype := :heatmap
        seriescolor --> :RdBu
        xcoords, ycoords, μ'
    end

    # Standard deviation section
    @series begin
        subplot := 2
        seriestype := :heatmap
        xcoords, ycoords, σ'
    end

    # Sources
    if !isnothing(sources)
        for subplot in (1, 2)
            @series begin
                seriestype := :scatter
                primary := false
                subplot := subplot
                markersize --> 3
                markercolor --> :orange
                markerstrokecolor --> :black
                markerstrokewidth -> 0.5
                label --> ""
                ((slice_name === :x ? sources.y : sources.x),
                    (slice_name === :z ? sources.y : sources.z))
            end
        end
    end

    # Receivers
    if !isnothing(receivers)
        for subplot in (1, 2)
            @series begin
                seriestype := :scatter
                primary := false
                subplot := subplot
                markersize --> 5
                markercolor --> :black
                markershape --> :dtriangle
                markerstrokecolor --> :white
                markerstrokewidth -> 0.5
                label --> ""
                ((slice_name === :x ? receivers.y : receivers.x),
                    (slice_name === :z ? receivers.y : receivers.z))
            end
        end
    end

    # Fix up annotations, aspect ratio, etc.
    # This is needed here because the scatter plots if used reset all
    # these.  It's cleaner therefore to just set them for an empty
    # series at the end.
    @series begin
        subplot := 1
        aspect_ratio := 1
        xguide --> xlabel
        yguide --> ylabel
        primary = false
        label := ""
        title --> "Mean $property_name"
        colorbar_title --> units
        [], []
    end

    @series begin
        subplot := 2
        aspect_ratio := 1
        xguide --> xlabel
        yguide --> ylabel
        clims --> (0, Inf)
        primary := false
        label := ""
        title --> "Standard deviation"
        colorbar_title --> units
        [], []
    end
end

end # module
