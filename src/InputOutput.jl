"""
# `IO`

Provides input and output for MCTomo files.
"""
module InputOutput

import CSV
import Tables

using StaticArrays: SVector

using ..Grids
using ..MCTomoTools: Chain, Model, Nodes, RawSample, RAW_SAMPLE_LEN_BYTES
using ..NameLists: read_namelist
using ..Settings: MCTomoSettings, getpath

export
    read_chains_burnin,
    read_initial_sample,
    read_receivers,
    read_raw_samples,
    read_settings,
    read_sources,
    read_temperatures,
    read_velocities!,
    read_velocities

"""
    Chain(settings::MCTomoSetting, chain_index) -> chain

Construct a chain by reading the raw samples from disk for chain
with index `chain_index`.

# Example
```
julia> settings = read_settings("MCTomo.inp")
MCTomoSettings(Dict{String, Dict{String, Dict{String,Any}}}(…), "sample_directory")

julia> chain_5 = Chain(settings, 5)
Chain{Float64,Vector{MCTomoTools.RawSample}}:
Model{Float64}:
 nodes: (184 nodes)
 body_noise0: [0.0, 0.0]
 body_noise1: [0.030349066418211227, 0.45723031209950116]
 surfacewave_noise0: [0.0214651962533745]
 surfacewave_noise1: [0.03637368954633682]
 locations: (61 sources)
Vector{MCTomoTools.RawSample}:
 steps: 2308000
```
"""
Chain(settings::MCTomoSettings, chain_index::Integer) =
    Chain(Model(settings, chain_index), read_raw_samples(settings, chain_index))

"""
    Grid([T=Float64,] settings::MCTomoSettings) -> grid

Create a new, undefined grid with coordinates and spacing
defined by `settings`, which holds the settings for an MCTomo run.
"""
function Grids.Grid(T, settings::MCTomoSettings)
    x0 = settings["grid_settings", "grid", "xmin"]
    x1 = settings["grid_settings", "grid", "xmax"]
    nx = settings["grid_settings", "grid", "nx"]
    y0 = settings["grid_settings", "grid", "ymin"]
    y1 = settings["grid_settings", "grid", "ymax"]
    ny = settings["grid_settings", "grid", "ny"]
    z0 = settings["grid_settings", "grid", "zmin"]
    z1 = settings["grid_settings", "grid", "zmax"]
    nz = settings["grid_settings", "grid", "nz"]
    Grids.Grid(T, x0, x1, nx, y0, y1, ny, z0, z1, nz)
end
Grids.Grid(settings::MCTomoSettings) = Grids.Grid(Float64, settings)

"""
    read_chains_burnin(io_or_file) -> chains::Dict{Int,Int}
    read_chains_burnin(settings::Settings.MCTomoSettings) -> chains::Dict{Int,Int}

Read the `chains.inp` file (either passed in as an `IO` object, file name,
or derived from the `settings` of a run), and return a dictionary `chains`
mapping the chain index to the number of burnin samples for each chain.
When sampling the chains, the burnin samples should be ignored.
"""
function read_chains_burnin(io::IO)
    nchains = parse(Int, readline(io))
    chains_burnin = Dict{Int,Int}()

    for linenumber in 2:(nchains+1)
        line = readline(io)
        tokens = split(line)
        length(tokens) == 2 || throw(ArgumentError(
            "unexpected number of columns ($(length(tokens)), not 2) " *
            "at line $linenumber"))
        chainid = parse(Int, tokens[1])
        burnin = parse(Int, tokens[2])
        haskey(chains_burnin, chainid) &&
            throw(ArgumentError("chain ID $chainid appears more than once"))
        chains_burnin[chainid] = burnin
    end

    chains_burnin
end
read_chains_burnin(file) = open(read_chains_burnin, file)
read_chains_burnin(settings::MCTomoSettings) =
    open(read_chains_burnin, joinpath(settings, "chains.inp"))

"""
    read_initial_sample(io_or_file) -> nodes::Nodes
    read_initial_sample(settings::MCTomoSettings, chain_index) -> nodes::Nodes

Read the initial sample node locations and properties
from `io_or_file`, returning a `Nodes` object.
"""
function read_initial_sample(io::IO)
    ncells = Int(read(io, Csize_t))
    coords = Vector{SVector{3,Cdouble}}(undef, ncells)
    properties = copy(coords)
    nodes = Nodes(coords, properties)
    for i in eachindex(nodes.coords)
        nodes.coords[i] = read(io, SVector{3,Cdouble})
        nodes.properties[i] = read(io, SVector{3,Cdouble})
    end
    nodes
end
read_initial_sample(file) = open(read_initial_sample, file)
read_initial_sample(settings::MCTomoSettings, chain_index) =
    read_initial_sample(joinpath(settings, "Results", "ini", "InitialSample_$(chain_index).dat"))

"""
    read_initial_sigma(io_or_file) -> sigma0_body, sigma1_body, sigma0_surface, sigma1_surface

Read the initial noise parameters for the body- and surface-wave
data from `io_or_file`, returning the noise parameters
as `Vector{Cdouble}`s.

!!! note
    The `InitialSigma_n.dat` files are plain-text files, not
    stream-accessed binary files.
"""
function read_initial_sigma(io::IO)
    line = readline(io)
    tokens = split(line)
    length(tokens) == 2 ||
        throw(ArgumentError("expected two columns on line 1 of initial sigma file"))
    nsigma_body = parse(Int, tokens[1])
    nsigma_surface = parse(Int, tokens[2])

    sigma0_body = Vector{Cdouble}(undef, nsigma_body)
    sigma1_body = copy(sigma0_body)
    sigma0_surface = Vector{Cdouble}(undef, nsigma_surface)
    sigma1_surface = copy(sigma0_surface)

    for i in eachindex(sigma0_body, sigma1_body)
        tokens = split(readline(io))
        sigma0_body[i] = parse(Cdouble, tokens[1])
        sigma1_body[i] = parse(Cdouble, tokens[2])
    end
    for i in eachindex(sigma0_surface, sigma1_surface)
        tokens = split(readline(io))
        sigma0_surface[i] = parse(Cdouble, tokens[1])
        sigma1_surface[i] = parse(Cdouble, tokens[2])
    end

    (; sigma0_body, sigma1_body, sigma0_surface, sigma1_surface)
end
read_initial_sigma(file) = open(read_initial_sigma, file)

"""
    read_receivers(io_or_file, T=Float64) -> x, y, z

Read the locations of the receivers from the file specified in `MCTomo.inp`,
usually called `breceivers.dat`, and return three `Vector{T}`s giving their
`x`, `y` and `z` coordinates in km.

Optionally specify the element type `T`.

---

    read_receivers(settings::MCTomoSettings, T=Float64) -> x, y, z

Read from the file specified within the `settings`.
"""
function read_receivers(io::IO, T=Float64)
    data = CSV.read(io, Tables.matrix; skipto=2, header=false, delim= ' ',
        ignorerepeated=true, types=T, strict=true)

    size(data, 2) == 3 || throw(
        ArgumentError("unexpected number of columns ($(size(data, 2)) not 3)"))

    x = map(T, @view(data[:,1]))
    y = map(T, @view(data[:,2]))
    z = map(T, @view(data[:,3]))
    (; x, y, z)
end
function read_receivers(settings::MCTomoSettings, T=Float64)
    read_receivers(getpath(settings, "likelihood_settings", "like_set", "breceivers_file"))
end
read_receivers(file, T=Float64) = open(io -> read_receivers(io, T), file)

"""
    read_raw_samples(io_or_file; n=nothing) -> raw_samples

Read `n` samples from a binary file, returning a `Vector{MCTomoTools.RawSample}`
`samples` containing:
- `step::Vector{Int}`: Sample number
- `accepted::Vector{Bool}`: Whether sample was accepted or not
- `vindex::Vector{Int}`: Index of vertex (or event?) adjusted
- `ncells::Vector{Int}`: Number of cells
- `misfit::Vector{T}`: Value of misfit
- `likelihood::Vector{T}`: Value of likelihood
- `x::Vector{T}`: x-coordinate of vertex (or event?)
- `y::Vector{T}`: y-coordinate of vertex
- `z::Vector{T}`: z-coordinate of vertex
- `vp::Vector{T}`: Vp value of vertex
- `vs::Vector{T}`: Vs value of vertex
- `density::Vector{T}`: density value of vertex
- `noise0::Vector{T}`: Linear noise parameter term (noise = noise0*distance + noise1)
- `noise1::Vector{T}`: Constant noise paramter term

By default (when `n` is `nothing`), all (remaining) samples in the stream
or file are read.

!!! note
    You should specify `n` if `io_or_file` is a stream which
    doesn't support `stat`, as this is used to determine the total number of
    samples in the file or `IO` object.
"""
function read_raw_samples(io::IO; n=nothing)
    samples_to_read = if isnothing(n)
        (stat(io).size - position(io))÷RAW_SAMPLE_LEN_BYTES
    else
        n
    end

    samples = [read(io, RawSample) for _ in 1:samples_to_read]
    if isnothing(n)
        !eof(io) && @warn("attempted to read all samples, but not at end-of-file")
    end

    samples
end
read_raw_samples(file; n=nothing) = open(io -> read_raw_samples(io; n), file)

"""
    read_raw_samples(settings::MCTomoSettings, chain_index; n=nothing) -> raw_samples

Read samples from chain `chain_index`, with the file location determined
by the path to the `MCTomo.inp` file in `settings`.
"""
function read_raw_samples(settings::MCTomoSettings, chain_index; n=nothing)
    file = joinpath(settings, "Results", "samples_$(chain_index).out")
    read_raw_samples(file; n=n)
end

"""
    read_sources(io_or_file, T=Float64) -> x, y, z, t

Read the locations of the sources from the file specified in `MCTomo.inp`,
usually called `bsources.dat`, and return four `Vector{T}`s giving their
`x`, `y` and `z` coordinates in km and time `t` in s.

Optionally specify the element type `T`.

---

    read_sources(settings::MCTomoSettings, T=Float64) -> x, y, z, t

Read from the file specified within the `settings`.
"""
function read_sources(io::IO, T=Float64)
    data = CSV.read(io, Tables.matrix; skipto=2, header=false, delim= ' ',
        ignorerepeated=true, types=T, strict=true)

    size(data, 2) == 4 || throw(
        ArgumentError("unexpected number of columns ($(size(data, 2)) not 4)"))

    x = map(T, @view(data[:,1]))
    y = map(T, @view(data[:,2]))
    z = map(T, @view(data[:,3]))
    t = map(T, @view(data[:,4]))
    (; x, y, z, t)
end
function read_sources(settings::MCTomoSettings, T=Float64)
    read_sources(getpath(settings, "likelihood_settings", "like_set", "bsources_file"))
end
read_sources(file, T=Float64) = open(io -> read_sources(io, T), file)

"""
    read_temperatures(io_or_file, T=Float64) -> temps
    read_temperatures(settings::MCTomoSettings, chain_index, T=Float64) -> temps

For a run with parallel tempering, read a file which describes how
individuals samples are swapped between chains and what their
temperatures are.

If using an `MCTomoSettings` instance, supply the index of the chain of
interest, `chain_index`.

`temps` is a named tuple of:
- `swapped::Vector{Bool}`: Whether each sample swapped chain
- `chain1_temp::Vector{T}`: Temperature of the sample's old chain, if swapped
- `chain2_temp::Vector{T}`: Temperature of the sample's new chain, if swapped
- `sample_temp::Vector{T}`: Temperature of the sample after swapping
"""
function read_temperatures(io::IO, T=Float64)
    data = CSV.read(io, Tables.matrix, header=false, delim=' ',
        ignorerepeated=true, types=T, strict=true)
    size(data, 2) == 4 || throw(
        ArgumentError("unexpected number of columns ($(size(data, 2)) not 4)"))
    swapped = map(Bool, @view(data[:,1]))
    chain1_temp = map(T, @view(data[:,2]))
    chain2_temp = map(T, @view(data[:,3]))
    sample_temp = map(T, @view(data[:,4]))
    (; swapped, chain1_temp, chain2_temp, sample_temp)
end
read_temperatures(file, T=Float64) = open(io -> read_temperatures(io, T), file)
read_temperatures(settings::MCTomoSettings, chain_index::Integer, T=Float64) =
    read_temperatures(joinpath(settings, "Results", "temperatures_$(chain_index).dat"))

"""
    read_velocities!(io_or_file, vp_grid, vs_grid) -> (vp_grid, vs_grid)

Read two velocity grids for Vp and Vs from `io_or_file`, which may
be an `IO` object (such as returned from `open` or the name of a file.
The `Grid`s `vp_grid` and `vs_grid` should have their dimensions
already defined, and their contents are replaced from the file.
"""
function read_velocities!(io::IO, vp_grid::Grid{T}, vs_grid::Grid{T}) where T
    size(vp_grid) == size(vs_grid) ||
        throw(ArgumentError("vp and vs grids must be the same size"))
    coordinates(vp_grid) == coordinates(vs_grid) ||
        throw(ArgumentError("vp and vs grids must have the same coordinates"))

    nx, ny, nz = size(vp_grid)

    data = CSV.read(io, Tables.matrix; header=false, delim=' ',
        ignorerepeated=true, types=T, strict=true)
    nlines, nz_file = size(data)
    
    iseven(nlines) || error("data do not contain an even number of lines")
    nz_file == nz || throw(ArgumentError("nz in file ($nx_file) ≠ nz in grid ($nz)"))
    npts = nlines÷2
    npts == nx*ny ||
        throw(ArgumentError("lateral number of points in file $npts ≠ $nx × $ny"))
    
    vp_data = reshape(data[1:end÷2,:], nx, ny, nz)
    vs_data = reshape(data[end÷2+1:end,:], nx, ny, nz)
    
    # Conversion happens here
    vp_grid.data .= vp_data
    vs_grid.data .= vs_data

    vp_grid, vs_grid
end

function read_velocities!(file, vp_grid::Grid{T}, vs_grid::Grid{T}) where T
    open(file) do io
        read_velocities!(io, vp_grid, vs_grid)
    end
end

"""
    read_velocities(io_or_file, settings::MCTomoSettings, T=Float64) -> vp::Grid, vs::Grid

Read two grids from `io_or_file` with the settings defined by the
`settings` read from an `"MCTomo.inp"` file.
"""
function read_velocities(io_or_file, settings::MCTomoSettings, T=Float64)
    x0 = settings["grid_settings", "grid", "xmin"]
    x1 = settings["grid_settings", "grid", "xmax"]
    nx = settings["grid_settings", "grid", "nx"]
    y0 = settings["grid_settings", "grid", "ymin"]
    y1 = settings["grid_settings", "grid", "ymax"]
    ny = settings["grid_settings", "grid", "ny"]
    z0 = settings["grid_settings", "grid", "zmin"]
    z1 = settings["grid_settings", "grid", "zmax"]
    nz = settings["grid_settings", "grid", "nz"]
    read_velocities(io_or_file, x0, x1, nx, y0, y1, ny, z0, z1, nz, T)
end

"""
    read_velocities(io_or_file, x0, x1, nx, y0, y1, ny, z0, z1, nz, T=Float64) -> vp::Grid, vs::Grid

Read two grids from `io_or_file`, given the coordinates defined
by the starting (`_0`) and ending (`_1`) coordinates and the number
of points (`n_`).
"""
function read_velocities(io, x0, x1, nx, y0, y1, ny, z0, z1, nz, T=Float64)
    vp = Grid(T, x0, x1, nx, y0, y1, ny, z0, z1, nz)
    vs = Grid(T, x0, x1, nx, y0, y1, ny, z0, z1, nz)
    read_velocities!(io, vp, vs)
end

"""
    read_velocities!(files, vp::Grid, vs::Grid) -> vps::Vector{Grid}, vss:Vector{Grid}

Read a set of `files` into vectors of `Grid`s, returning two vectors:
one for Vp and one for Vs.
"""
function read_velocities!(files::AbstractArray{<:AbstractString},
        vp::Grid{T}, vs::Grid{T}; verbose=false) where T
    vps = Grid{T}[]
    vss = similar(vps)
    for file in files
        verbose && print(stderr, "\rReading file $(basename(file))")
        read_velocities!(file, vp, vs)
        push!(vps, copy(vp))
        push!(vss, copy(vs))
    end
    verbose && println(stderr)
    vps, vss
end

"""
    read_settings(io_or_file) -> ::Settings.MCTomoSettings

Read the `MCTomo.inp` input file and return a [`Settings.MCTomoSettings`](@ref)
structure.

If `io_or_file` is an `IO` object, then it is assumed that the
MCTomo run files are in the current directory.
"""
function read_settings(io::IO; file=nothing)
    dir = isnothing(file) ? "." : dirname(file)
    MCTomoSettings(read_namelist(io; file=file), dir)
end
read_settings(file) = open(io -> read_settings(io; file=file), file)

end # module
