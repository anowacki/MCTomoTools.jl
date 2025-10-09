function MCTomoTools.plot_section_makie(
    chains::Union{Chain,AbstractArray{<:Chain}},
    grid::Grid,
    property::Symbol;
    sources=nothing,
    receivers=nothing,
    burnin=0,
    thin=1000,
    x=nothing, y=nothing, z=nothing,
    method=:brute,
    figure=(),
    axis=(),
    axis_stdev=axis,
    heatmap=(),
    heatmap_stdev=(),
)
    count(isnothing, (x, y, z)) == 2 ||
        throw(ArgumentError("one and only one of x, y or z must be given"))

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
    μ, σ = MCTomoTools.sample_average_grid(chains, slice, property; burnin, thin, method)
    μ = dropdims(μ, dims=slice_dim)
    σ = dropdims(σ, dims=slice_dim)

    # Work out correct axes and things for this slice
    xcoords = (slice_name === :x ? grid.y : grid.x)
    ycoords = (slice_name === :z ? grid.y : grid.z)
    xlabel = (slice_name === :x ? "Northing / km" : "Easting / km")
    ylabel = (slice_name === :z ? "Northing / km" : "Depth / km")

    property_name = uppercasefirst(String(property))
    units = property === :density ? "g/cm^3" : "km/s"

    fig = Makie.Figure(; figure...)

    # Mean plot
    ax_mean, hm_mean = Makie.heatmap(
        fig[1,1], xcoords, ycoords, μ;
        axis=(;
            xlabel, ylabel, aspect=Makie.DataAspect(),
            title="Mean $property_name",
            axis...
        ),
        colormap=:RdBu,
        heatmap...
    )
    cbar_mean = Makie.Colorbar(fig[1,2], hm_mean; label="$property_name ($units)")

    # Source and receiver locations if any
    _plot_sources!(ax_mean, sources, slice_name)
    _plot_receivers!(ax_mean, receivers, slice_name)

    # Standard deviation plot
    ax_stdev, hm_stdev = Makie.heatmap(
        fig[2,1], xcoords, ycoords, σ;
        axis=(
            xlabel, ylabel, aspect=Makie.DataAspect(),
            title="Standard deviation ($property_name)",
            axis_stdev...
        ),
        colormap=:grays,
        heatmap_stdev...
    )
    cbar_stdev = Makie.Colorbar(fig[2,2], hm_stdev; label="Stdev $property_name ($units)")
    _plot_sources!(ax_stdev, sources, slice_name)
    _plot_receivers!(ax_stdev, receivers, slice_name)

    fig
end

function _plot_sources!(ax, sources, slice_name; kwargs...)
    if !isnothing(sources)
        Makie.scatter!(
            ax,
            (slice_name === :x ? sources.y : sources.x),
            (slice_name === :z ? sources.y : sources.z);
            color=:orange,
            marker=:circle,
            strokecolor=:white,
            strokewidth=1,
            kwargs...
        )
    end
end

function _plot_receivers!(ax, receivers, slice_name; kwargs...)
    if !isnothing(receivers)
        Makie.scatter!(
            ax,
            (slice_name === :x ? receivers.y : receivers.x),
            (slice_name === :z ? receivers.y : receivers.z);
            marker=:dtriangle,
            strokecolor=:white,
            strokewidth=1,
            color=:darkgreen,
        )
    end
end
