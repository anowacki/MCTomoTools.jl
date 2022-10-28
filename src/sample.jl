# Create samples from a chain

"""
    sample_chain(func, chain::Chain; burnin=0, thin=0) -> nsamples

For every `thin`th model in a chain of proposals `steps`, run `func`,
optionally discarding the first `burnin` samples.

`func(sample)` is a function taking a single argument, `sample`, which is
a named tuple `(model, isample)`. `model` is the `Model` at each
(thinned) step, while `isample` is the index of the sample, accounting
for burn-in and thinning.  Note that the same `model` passed in to this function
is updated in place and then passed back to `func`, so you may need to
`copy(model)` within `func` if you do not wish to store aliased copies
of `model`.

# Examples
Calculate the mean number of nodes for all models in the chain with no burn-in.
```
julia> mean_nnodes = let sum_nnodes = 0
           nsamples = sample_chain!(model, steps; thin=1) do sample
               sum_nnodes += length(sample.model.nodes)
           end
           sum_nnodes/nsamples
       end
```
"""
sample_chain(func, chain::Chain; kwargs...) =
    sample_chain!(func, initial_model(chain), chain.samples; kwargs...)

"""
    sample_chain(func, model::Model, steps::Vector{RawSample}; burnin=0, thin=0) -> nsamples

As above, but using an individual `Model` and set of `steps`.
"""
sample_chain(func, model::Model, steps::AbstractVector{<:RawSample}; kwargs...) =
    sample_chain!(func, deepcopy(model), steps; kwargs...)

"""
    sample_chain!(func, model::Model, steps::Vector{RawSample}; burnin=0, thin=1000) -> nsamples

As for [`sample_chain`](@ref), but **modifies `model` in place**.  If
you want to avoid this, copy `model` before passing it in.
"""
function sample_chain!(func, model::Model, steps::AbstractVector{<:RawSample};
    burnin=0, thin=1000
)
    nmodels = 0
    for (istep, step) in pairs(steps)
        update!(model, step, istep)
        if _add_thinned_step(istep, burnin, thin)
            nmodels += 1
            sample = (isample=nmodels, model=model)
            func(sample)
        end
    end
    nmodels
end

"""
    sample_grid!(grid::Grid, chain::Chain, property::Symbol; kwargs...) -> grids::Vector{Grid}

For a `chain` of steps, perturb the initial model by all the accepted proposals
steps, and compute the grid corresponding to `property` at every
`thin`th accepted step.  Optionally discard the first `burnin` accepted
proposals.  A `Vector{Grid{T}}` is returned.  `grid` is used as a buffer and
modified.

`property` may be one of:
- `:vp`: P-wave velocity (km/s)
- `:vs`: S-wave velocity (km/s)
- `:density`: Density (g/cm³) (but note density is not independently sampled
currently in MCTomo)

`grid` is a [`Grid`](@ref) of the desired size and coordinates which will be
filled from the voronoi nodes in `model`.  This can be constructed to be the
same as the travel time computation grid used in the inversion by calling
[`settings = read_settings("MCTomo.inp")`](@ref read_settings) then
`grid = Grid(settings)`.  However, you may want to use a coarser or smaller
grid for the final set of models to speed up resampling from the saved chain.

!!! warn
    If `thin` is too small, there is a chance that all of your memory will be
    used up by constructing so many `Grid`s!
"""
sample_grid!(grid::Grid, chain::Chain, property::Symbol; kwargs...) =
    sample_grid!(grid, initial_model(chain), chain.samples, property; kwargs...)

"""
    sample_grid!(grid::Grid{T}, model::Model, steps::Vector{RawSample}, property::Symbol; burnin=0, thin=1000, method=:brute)

As above, but using an individual `model` and `steps`.  Note that **both grid
and model are updated in place** by this method.  Use [`sample_grid`](@ref)
if that is not desired.

`steps` may be read by calling [`read_raw_samples`](@ref).
"""
function sample_grid!(grid::Grid, model::Model,
    steps::AbstractVector{<:RawSample}, property::Symbol; burnin=0, thin=1000, method=:brute
)
    property_index = _property_index(property)

    grids = typeof(grid)[]

    for (istep, step) in pairs(steps)
        update!(model, step, istep)
        if istep > burnin && (istep - burnin)%thin == 0
            Neighbours.fill_from_nodes!(grid, model.nodes, property_index; method)
            push!(grids, deepcopy(grid))
        end
    end
     
    grids
end

"""
    sample_grid!(func, grid::Grid, chain::Chain, property; burnin=0, thin=1000, method=:brute) -> nmodels::Int

At each model in the chain, run `func`, and return the number of models
evaluated.

`func` should be a function which takes two input arguments and
returns `nothing`.  The first argument is the current `model::Model` at
each step in the chain, and the second argument is the `grid::Grids.Grid`
filled in from the model at that step.  `grid` is a reference to the
original grid passed in, so you should call `copy(grid)` if you don't
want to simply alias the original input.

`property` may be one of:
- `:vp`: P-wave velocity (km/s)
- `:vs`: S-wave velocity (km/s)
- `:density`: Density (g/cm³) (but note density is not independently sampled
currently in MCTomo)

Optionally pass `:kdtree` as `method` to use a KD tree to perform the nearest-neighbour
search in constructing the grid.  The default (`:brute`) performs a brute-force
search which may be quicker for small grids or 2D slices, but may be slower for
full 3D grids.

# Examples
To copy the functionality of the first form of `sample_grid!`, you could
do:
```
julia> grids = typeof(grid)[];

julia> sample_grid!(model, grid, steps, property) do model, grid
           push!(grids, copy(grid))
       end
```

To calculate the 1D mean profile, you could do:
```
julia> using Statistics

julia> profile = zeros(eltype(grid), axes(grid, 3));

julia> nsamples = sample_grid!(model, grid, steps, property) do model, grid
           for ilayer in eachindex(profile, grid.z)
               profile[ilayer] += mean(@view(grid[:,:,ilayer]))
           end
       end

julia> profile ./= nsamples
```
"""
sample_grid!(func, grid::Grid, chain::Chain, property::Symbol; kwargs...) =
    sample_grid!(func, grid, initial_model(chain), chain.samples, property; kwargs...)

"""
    sample_grid!(func, grid, model, steps, property; burnin=0, thin=1000, method=:brute) -> nsamples

As above, but using individual an `model` and `steps`.  **`model` is
modified in place** by this method; use [`sample_grid`](@ref) if that
is not wanted.
"""
function sample_grid!(func::F, grid::Grid, model::Model, steps::AbstractVector{<:RawSample},
    property::Symbol; burnin=0, thin=1000, method=:brute
) where F
    property_index = _property_index(property)

    nmodels = 0
    for (istep, step) in pairs(steps)
        update!(model, step, istep)
        if istep > burnin && (istep - burnin)%thin == 0
            nmodels += 1
            Neighbours.fill_from_nodes!(grid, model.nodes, property_index; method)
            func(model, grid)
        end
    end
    nmodels
end

"""
    sample_grid(func, chain::Chain, grid::Grid, property; burnin=0, thin=1000, method=:brute) -> nsamples
    sample_grid(chain::Chain, grid::Grid, property; burnin=0, thin=1000, method=:brute) -> ::Vector{Grids.Grid}

As for [`sample_grid!`](@ref), but without modifying `grid`.
"""
sample_grid(func, chain::Chain, grid::Grid, property::Symbol; kwargs...) =
    sample_grid!(func, deepcopy(grid), initial_model(chain), chain.samples, property; kwargs...)
sample_grid(chain::Chain, grid::Grid, property::Symbol; kwargs...) =
    sample_grid!(deepcopy(grid), initial_model(chain), chain.samples, property; kwargs...)

"""
    sample_grid(func, model::Model, steps::Vector{RawSample}, grid, property; burnin=0, thin=1000, method=:brute)
    sample_grid(model::Model, steps::Vector{RawSample}, grid, property; burnin=0, thin=1000, method=:brute) -> grids::Vector{Grids.Grid}

As above, but using an individual `model` and `steps`
"""
sample_grid(model, steps, grid, property; kwargs...) =
    sample_grid!(deepcopy(model), steps, deepcopy(grid), property; kwargs...)
sample_grid(func, model, steps, grid, property; kwargs...) = 
    sample_grid!(func, deepcopy(model), steps, deepcopy(grid), property; kwargs...)

"""
    sample_average_grid!(grid::Grid, chain::Chain, property::Symbol; burnin=0, thin=1000, slowness=false, method=:brute) -> μ::Grid, σ::Grid

As for [`sample_average_grid!`](@ref), except `grid` is modified during
the samping process and used as a buffer.
"""
sample_average_grid!(grid::Grid, chain::Chain, property::Symbol; kwargs...) =
    sample_average_grid!(grid, initial_model(chain), chain.samples, property; kwargs...)

"""
    sample_average_grid!(grid::Grid, model::Model, steps::Vector{RawSample}, property; burnin=0, thin=1000, method=:brute) -> μ::Grid, σ::Grid

As above, but using an individual `Model` and `steps`.  Note that **`model`
is modified in place** in this method.
"""
function sample_average_grid!(grid::Grid, model::Model,
    steps::AbstractArray{<:RawSample}, property::Symbol;
    burnin=0, thin=1000, slowness=false, method=:brute
)
    property == :density && throw(ArgumentError("density averaging not implemented"))
    property_index = _property_index(property)
    
    # Note: these are `Array`s not `Grid`s as we have not
    #       written `Base.similar` or `Base.zero` methods for `Grid`,
    #       but that is what we want.
    # Running total velocity/slowness and eventually mean
    # Running squared velocity/slowness and eventually variance
    μ = zero(grid)
    σ² = zero(grid)
    # Temporary grid for running calculations
    temp = similar(grid)

    nmodels = 0
    for (istep, step) in pairs(steps)
        update!(model, step, istep)
        if istep > burnin && (istep - burnin)%thin == 0

            # Average slowness if requested
            nodes = slowness ? _to_slowness(model.nodes) : model.nodes
            Neighbours.fill_from_nodes!(grid, nodes, property_index; method)

            nmodels = _running_mean_stdev!(temp, nmodels, μ, σ², grid)
        end
    end

    μ, σ = _finalise_mean_stdev!(nmodels, μ, σ²)

    # Convert to velocity if needed and put back in `Grid`s
    slowness ?
        (Grid(grid, inv.(μ)), Grid(grid, inv.(σ))) :
        (Grid(grid, μ), Grid(grid, σ))
end

"""
    sample_average_grid(chain, grid, property; burnin=0, thin=1000, method=:brute) -> μ::Grid, σ::Grid

For a `chain` of samples, perturb the inital model by all the accepted
proposal steps and compute the mean (`μ`) and standard deviation (`σ`) of the
grids corresponding to `property` at every `thin`th accepted step.
Optionally discard the first `burnin` accepted proposals.  If `slowness`
is `true`, then calculate the harmonic mean (inverse of the mean slownesses)
instead.

`steps` may be read by calling [`read_raw_samples`](@ref).

`property` may be one of:
- `:vp`: P-wave velocity (km/s)
- `:vs`: S-wave velocity (km/s)
- `:density`: Density (g/cm³) (but note density is not independently sampled
currently in MCTomo)

`grid` is a [`Grid`](@ref) of the desired size and coordinates, a copy of
which will be filled from the voronoi nodes in `model`.  This can be constructed
to be the same as the travel time computation grid used in the inversion by calling
[`settings = read_settings("MCTomo.inp")`](@ref read_settings) then
`grid = Grid(settings)`.  However, you may want to use a coarser or smaller
grid for the final set of models to speed up resampling from the saved chain.

Optionally pass `:kdtree` as `method` to use a KD tree to perform the nearest-neighbour
search in constructing the grid.  The default (`:brute`) performs a brute-force
search which may be quicker for small grids or 2D slices, but may be slower for
full 3D grids.
"""
sample_average_grid(chain::Chain, grid::Grid, property::Symbol; kwargs...) =
    sample_average_grid!(deepcopy(grid), initial_model(chain), chain.samples, property; kwargs...)

"""
    sample_average_grid(model, steps, grid, property; burnin=0, thin=1000. method=:brute) -> μ::Grid, σ::Grid

As above, but directly pass a `Model::model` and set of `steps` as a
`Vector{RawSample}`.
"""
sample_average_grid(model::Model, steps::AbstractArray{<:RawSample}, grid::Grid, property::Symbol; kwargs...) =
    sample_average_grid!(deepcopy(grid), deepcopy(model), steps, property; kwargs...)

"""
    sample_average_grid(chains::Vector{Chain}, grid::Grid, property::Symbol; burnin=0, thin=1000, slowness=false, method=:brute) -> μ::Grid, σ::Grid

Calculate the combined mean `μ` and standard deviation `σ` for several `chains`.
"""
function sample_average_grid(
    chains::AbstractArray{<:Chain}, grid::Grid, property::Symbol;
    burnin=0, thin=1000, slowness=false, method=:brute
)
    # Temporary grid for running calculations
    buffer_grid = Grid(grid, similar(grid))

    # Note: these are `Array`s not `Grid`s as we have not
    #       written `Base.similar` or `Base.zero` methods for `Grid`,
    #       but that is what we want.
    # Mean grid for each chain
    μ = zero(grid)
    # Standard deviation for each chain
    σ = zero(grid)

    total_nsamples = 0
    for (ichain, chain) in pairs(chains)
        model = initial_model(chain)
        steps = chain.samples
        nsamples = length(burnin:thin:length(chain))
        μᵢ, σᵢ = sample_average_grid!(buffer_grid, model, steps, property;
            burnin=burnin, thin=thin, slowness=slowness, method=method)
        
        # Compute new σ in place first because we need both the old and new means
        σ .= _combined_stdev.(μ, σ, total_nsamples, μᵢ, σᵢ, nsamples)
        # Can now fill in the old means with the combined means
        μ .= _combined_mean.(μ, total_nsamples, μᵢ, nsamples)
    end

    Grid(grid, μ), Grid(grid, σ)
end

"Calculate the running mean and standard devation using Welford's algorithm"
function _running_mean_stdev!(buf, count, μ, σ², new)
    count += 1
    delta = buf
    delta .= new .- μ
    μ .+= delta./count
    σ² .+= delta.*(new .- μ)
    count
end

"""Compute the final mean and standard deviation from the running totals
returned by [`_running_mean_stdev!`](@ref MCTomoTools._running_mean_stdev!)."""
function _finalise_mean_stdev!(count, μ, σ²)
    σ² .= sqrt.(σ²./count)
    σ = σ²
    μ, σ
end

"""Standard deviation of two samples, each with mean μᵢ, standard deviation
σᵢ and sample size nᵢ."""
function _combined_stdev(μ₁, σ₁, n₁, μ₂, σ₂, n₂)
    sqrt((n₁*σ₁ + n₂*σ₂)/(n₁ + n₂) + n₁*n₂*(μ₁ - μ₂)^2/(n₁ + n₂))
end

"Mean of two samples, each with mean μᵢ and sample size nᵢ."
_combined_mean(μ₁, n₁, μ₂, n₂) = (n₁*μ₁ + n₂*μ₂)/(n₁ + n₂)


"""
    _to_slowness(old::Nodes) -> new::Nodes

Fill the nodes in `new` with the slownesses and inverse densities of `old`.
"""
_to_slowness(old::Nodes) = Nodes(old.coords, [1 ./ p for p in old.properties])

"""
    sample_source!(model::Model, steps::Vector{RawSample}, source_index; burnin=0, thin=1000) -> locations::Vector

For a starting `model`, perturb the model by all the accepted proposals steps
in `steps`, and compute the coordinates corresponding to one of the body-wave
sources at every `thin`th accepted step.  Optionally discard the first `burnin`
accepted proposals.  A `Vector` of locations is returned.

Note that each element of the returned vector is a named tuple of:
- `x`: Easting in km
- `y`: Northing in km
- `z`: Depth in km
- `t`: Difference from original origin time in s

`steps` may be read by calling [`read_raw_samples`](@ref).
"""
function sample_source!(model::Model, steps::AbstractArray{<:RawSample}, source_index::Integer;
    burnin=0, thin=1000
)
    first_coord = first(model.locations)
    coord = (x=first_coord[1], y=first_coord[2], z=first_coord[3], t=first_coord[4])
    coords = typeof(coord)[]

    for (istep, step) in pairs(steps)
        update!(model, step, istep)
        if istep > burnin && (istep - burnin)%thin == 0
            x, y, z, t = model.locations[source_index]
            push!(coords, (; x, y, z, t))
        end
    end

    coords
end

"""
    sample_source(model, steps, source_index; kwargs...)

As for [`sample_sources!`](@ref), but does not mutate its input arguments.
"""
sample_source(model, steps, source_index; kwargs...) =
    sample_source!(deepcopy(model), steps, source_index; kwargs...)

"""
    sample_sources!(model::Model, steps::AbstractArray{<:RawSample}; thin=1000)

For a starting `model`, perturb the model by all the acepted proposal steps
in `steps`, and compute the coordinates corresponding to all of the body-wave
sources at every `thin`th accepted step.  Optionally discard the first `burnin`
accepted samples.  A `Vector{Vector}}` of locations is returned, with the first
`Vector` corresponding to the first source, the second to the second, and so on.

Note that each element of each vector of the returned vector is a named tuple of:
- `x`: Easting in km
- `y`: Northing in km
- `z`: Depth in km
- `t`: Difference from original origin time in s

`steps` may be read by calling [`read_raw_samples`](@ref).
"""
function sample_sources!(model::Model, steps::AbstractArray{<:RawSample}; burnin=0, thin=1000)
    first_coord = first(model.locations)
    coord = (x=first_coord[1], y=first_coord[2], z=first_coord[3], t=first_coord[4])
    coords = [typeof(coord)[] for _ in eachindex(model.locations)]

    for (istep, step) in pairs(steps)
        update!(model, step, istep)
        if istep > burnin && (istep - burnin)%thin == 0
            for isource in eachindex(coords, model.locations)
                x, y, z, t = model.locations[isource]
                push!(coords[isource], (; x, y, z, t))
            end
        end
    end

    coords
end

"""
    sample_sources(model, steps; thin=1000)

As for [`sample_sources!`](@ref), but does not modify its input arguments.
"""
sample_sources(model, steps; thin=1000) = sample_sources!(deepcopy(model), steps; thin=thin)

"""
    _property_index(property::Symbol) -> index::Int

Convert `property` into an `index` suitable for indexing into
the 
"""
function _property_index(property::Symbol)
    if property === :vp
        1
    elseif property === :vs
        2
    elseif property === :density
        3
    else
        throw(ArgumentError("property must be `:vp`, `:vs` or `:density`"))
    end
end

"Helper function to choose whether to include a step when sampling from a chain"
@inline function _add_thinned_step(istep, burnin, thin)
    istep > burnin && (istep - burnin)%thin == 0
end
