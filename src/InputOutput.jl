"""
# `IO`

Provides input and output for MCTomo files.
"""
module InputOutput

import CSV
import Tables

using ..Grids

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

end # module
