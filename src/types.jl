# MCTomoTools types

"""
    Nodes{Ndim,Nprop,Tdim,Tprop}

Struct holding a set of nodes of dimension `Ndim` and element type `T`.
`Nprop` gives the number of properties held at each node, and `U` is the
element type of the properites.
Each node's coordinates are held as an `SVector{Tdim,Ndim}`, and properties
are an `SVector{Tprop,Nprop}`
"""
struct Nodes{Ndim,Nprop,Tdim,Tprop}
    coords::Vector{SVector{Ndim,Tdim}}
    properties::Vector{SVector{Nprop,Tprop}}
end

Base.size(nodes::Nodes) = size(nodes.coords)
Base.length(nodes::Nodes) = length(nodes.coords)
Base.getindex(nodes::Nodes, i::Int) = nodes.coords[i]
Base.setindex!(nodes::Nodes, v, i::Int) = nodes.coords[i] = v
Base.lastindex(nodes::Nodes) = lastindex(nodes.coords)
Base.axes(nodes::Nodes) = axes(nodes.coords)
Base.copy(nodes::Nodes) = Nodes(copy(nodes.coords), copy(nodes.properties))

function Base.push!(nodes::Nodes, coords, properties)
    push!(nodes.coords, coords)
    push!(nodes.properties, properties)
    nodes
end

function Base.deleteat!(nodes::Nodes, i::Int)
    deleteat!(nodes.coords, i)
    deleteat!(nodes.properties, i)
    nodes
end

"""
    Model

Instance of a model in MCTomo.

A model holds all the information necessary to describe a single state
of the model at a point in the Markov chain.  In MCTomo, models are
parameterised by a number of `nodes` which hold the position and
velocities (and density) of the Voronoi nodes of the space.  It also
holds the `locations` of each source.
"""
struct Model{T}
    nodes::Nodes{3,3,T,T}
    body_noise0::Vector{T}
    body_noise1::Vector{T}
    surfacewave_noise0::Vector{T}
    surfacewave_noise1::Vector{T}
    locations::Vector{SVector{4,T}}
end

"""
    Model(T) -> ::Model{T}
    Model() -> ::Model{$Cdouble}

Construct an empty `Model`, optionally with element type `T`.
"""
Model(T) = Model(Nodes(SVector{3,T}[], SVector{3,T}[]), T[], T[], T[], T[], SVector{4,T}[])
Model() = Model(Cdouble)

"""
    Model(settings::Settings.MCTomoSettings, chain_index) -> ::Model{$Cdouble}

Create an initial model from `settings`, which holds the settings for an
MCTomo run.  This is done for the single chain with index `chain_index`
"""
function Model(settings::Settings.MCTomoSettings, chain_index::Integer)
    ini_dir = joinpath(settings, "Results", "ini")
    noise = InputOutput.read_initial_sigma(joinpath(ini_dir, "InitialSigma_$(chain_index).dat"))
    nodes = InputOutput.read_initial_sample(joinpath(ini_dir, "InitialSample_$(chain_index).dat"))
    sources = InputOutput.read_sources(settings)
    Model(nodes, noise..., [@SVector[x, y, z, t] for (x, y, z, t) in zip(sources...)])
end

function Base.copyto!(new::Model, old::Model)
    if length(new.nodes) != length(old.nodes)
        resize!(new.nodes.coords, length(old.nodes))
        resize!(new.nodes.properties, length(old.nodes))
    end
    new.nodes.coords .= old.nodes.coords
    new.nodes.properties .= old.nodes.properties

    if length(new.body_noise0) != length(old.body_noise0)
        resize!(new.body_noise0)
        resize!(new.body_noise1)
    end
    new.body_noise0 .= old.body_noise0
    new.body_noise1 .= old.body_noise1

    if length(new.surfacewave_noise0) != length(old.surfacewave_noise0)
        resize!(new.surfacewave_noise0)
        resize!(new.surfacewave_noise1)
    end
    new.surfacewave_noise0 .= old.surfacewave_noise0
    new.surfacewave_noise1 .= old.surfacewave_noise1

    if length(new.locations) != length(old.locations)
        resize!(new.locations)
    end
    new.locations .= old.locations

    new
end

function Base.show(io::IO, ::MIME"text/plain", model::Model{T}) where T
    print(io, """Model{$T}:
         nodes: ($(length(model.nodes)) nodes)
         body_noise0: $(model.body_noise0)
         body_noise1: $(model.body_noise1)
         surfacewave_noise0: $(model.surfacewave_noise0)
         surfacewave_noise1: $(model.surfacewave_noise1)
         locations: ($(length(model.locations)) sources)
        """)
end

"""
    RawSample

Struct holding the raw information held in a binary MCTomo sample chain
file, typically called `samples_n.out`.

A sample may represent one of several different `step`s, indicated
by an integer which maps to the kind of update (per [`STEP_INDEX`](@ref)):

|`step`|Type                |`vindex`          |Updated fields              |
|:-----|:-------------------|:-----------------|:---------------------------|
|`1`   |Cell birth          |Always `0`        |`ncells,x,y,z,vp,vs,density`|
|`2`   |Cell death          |Old node index    |`ncells`                    |
|`3`   |Cell move           |Node index        |`x,y,z`                     |
|`4`   |Velocity change     |Node index        |`vp,vs,density`             |
|`5`   |Body noise change   |P (`1`) or S (`1`)|`noise0, noise1`            |
|`6`   |Surface noise change|?                 |`noise0, noise1`            |
|`7`   |Source move         |Source index      |`x,y,z,vp` (`vp` is time)   |

`RawSample`s can be read from disk or an `IO` object with a call to
`read(io, RawSample)`.
"""
struct RawSample
    step::Cint
    # Written as a Fortran iso_c_binding: c_bool, which is defined
    # as a C99 _Bool.  Usually this is 1 byte, but is system dependent!
    accepted::Int8
    vindex::Cint
    ncells::Csize_t
    misfit::Cdouble
    unweighted_misfit::Cdouble
    likelihood::Cdouble
    x::Cdouble
    y::Cdouble
    z::Cdouble
    vp::Cdouble
    vs::Cdouble
    density::Cdouble
    noise0::Cdouble
    noise1::Cdouble
end

@generated function Base.read(io::IO, ::Type{RawSample})
    names = fieldnames(RawSample)
    types = fieldtypes(RawSample)
    quote
        $([:($field = read(io, $typ)) for (field, typ) in zip(names, types)]...)
        RawSample($([:($f) for f in names]...))
    end
end

# Need to do this as otherwise `.accepted` is written as an Int32,
# since structs are not packed
@generated function Base.write(io::IO, sample::RawSample)
    names = fieldnames(RawSample)
    quote
        n = 0
        $([:(n += write(io, sample.$(field))) for field in names]...)
        n
    end
end

function Base.show(io::IO, ::MIME"text/plain", sample::RawSample)
    print(io, """RawSample:
         step: $(sample.step) ($(STEP_DESCRIPTION[sample.step]))
         accepted: $(Bool(sample.accepted))
         vindex: $(Int(sample.vindex))
         ncells: $(Int64(sample.ncells))
         misfit: $(sample.misfit)
         unweighted_misfit: $(sample.unweighted_misfit)
         likelihood: $(sample.likelihood)
         x: $(sample.x)
         y: $(sample.y)
         z: $(sample.z)
         vp: $(sample.vp)
         vs: $(sample.vs)
         density: $(sample.density)
         noise0: $(sample.noise0)
         noise1: $(sample.noise1)
        """)
end

"Mapping of the value of `step` in `RawSample` to the kind of proposed change"
const STEP_DESCRIPTION = Dict(
    1 => :cell_birth,
    2 => :cell_death,
    3 => :cell_move,
    4 => :velocity_change,
    5 => :bsigma_change,
    6 => :ssigma_change,
    7 => :loc_change,
)

"Mapping of the name of a proposed change to its `step` index in a `RawSample`"
const STEP_INDEX = Dict(desc => index for (index, desc) in STEP_DESCRIPTION)

"Length in bytes of the `RawSample` struct on disk"
const RAW_SAMPLE_LEN_BYTES = sum(sizeof, fieldtypes(RawSample))

"""
    Chain{T,V}

Struct holding both an initial model and all the steps within a
Markov chain which perturb it.
"""
struct Chain{T,V<:AbstractVector{<:RawSample}}
    initial_model::Model{T}
    _buffer_model::Model{T}
    samples::V
end

"""
    Chain(model::Model, samples::AbstractVector{RawSample}) -> chain

Construct a new `Chain` from a `model` and a set of `samples` (steps
in the Markov chain).

The following accessors are defined for `Chain`:
- `initial_model(chain)`: Returns the starting model of the chain, which
  chain then be modified by `MCTomoTools.update!` if necessary.

`Chain`s support the following methods which access the individual steps
in the chain, or return the length of the chain:
- `getindex(chain, i)`/`chain[i]`: Get the `i`th step
- `eachindex(chain)`: The indices of the chain
- `length(chain[[, burnin=0], thin=1])`: The length of a chain when
  sampled each `thin` steps after discarding the first `burnin` steps
"""
Chain(model::Model, samples::AbstractVector{<:RawSample}) =
    Chain(model, deepcopy(model), samples)

"""
    initial_model(chain::Chain) -> model::Model

Return the starting model of the `chain`, suitable for mutation.
"""
initial_model(chain::Chain) = copyto!(chain._buffer_model, chain.initial_model)

Base.getindex(chain::Chain, i::Int) = chain.samples[i]
Base.eachindex(chain::Chain) = eachindex(chain.samples)
Base.length(chain::Chain) = length(chain.samples)
Base.length(chain::Chain, thin::Integer) = length(chain, 0, thin)
Base.length(chain::Chain, burnin::Integer, thin::Integer) =
    (length(chain) - burnin)÷thin
Base.size(chain::Chain) = size(chain.samples)

function Base.show(io::IO, mime::MIME"text/plain", chain::Chain{T,V}) where {T,V}
    println(io, "Chain{$T,$V}:")
    show(io, mime, initial_model(chain))
    println(io, "$V:")
    print(io, " steps: $(length(chain))")
end
