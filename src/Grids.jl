"""
# `Grids`

Module for dealing with velocity grids
"""
module Grids

export
    Grid,
    MultiGrid,
    coordinates

"""
    Grid{T,R} <: AbstractArray{T,3}

Grid containing a 3D volume of points spaced equally.

# User-accessible fields
- `x`, `y` and `z`: Ranges containing the x, y and z coordinates in
  each dimension.  Units km.
- `grid`: 3D array containing the values at each point.  Units km/s.

# Indexing
`Grid`s are `AbstractArray`s, and so can be indexed just like normal
3D arrays with `getindex` and `setindex!` (e.g., `grid[:,:,3]`).
"""
struct Grid{T,R} <: AbstractArray{T,3}
    x::R
    y::R
    z::R
    data::Array{T,3}
end

"""
    Grid(x0, x1, y0, y1, z0, z1, grid) -> ::Grid

Construct a new `Grid` giving the start and end coordinates (`x0`
and `x1`) in each of the three dimensions, and a 3D array of values
`grid`.  The number of points in `grid` determines how many
coordinate points there are between e.g. `x0` and `x1`.
"""
function Grid(x0, x1, y0, y1, z0, z1, data::AbstractArray{T,3}) where T
    nx, ny, nz = size(grid)
    # Convert coordinates to correct precision
    x = map(T, range(x0, stop=x1, length=nx))
    y = map(T, range(y0, stop=y1, length=ny))
    z = map(T, range(z0, stop=z1, length=nz))
    Grid(x, y, z, data)
end

"""
    Grid(T, x0, x1, nx, y0, y1, ny, z0, z1, nz) -> ::Grid{T}
    Grid(x0, x1, nx, y0, y1, ny, z0, z1, nz) -> ::Grid{Float64}

Create a new empty `Grid` by specifying the start and end values
of the coordinates in each dimension, and the number of grid points.
Specify the element type `T`, which defaults to `Float64` if not given.
"""
function Grid(T, x0, x1, nx, y0, y1, ny, z0, z1, nz)
    x = map(T, range(x0, stop=x1, length=nx))
    y = map(T, range(y0, stop=y1, length=ny))
    z = map(T, range(z0, stop=z1, length=nz))
    data = Array{T}(undef, nx, ny, nz)    
    Grid(x, y, z, data)
end
Grid(x0, x1, nx, y0, y1, ny, z0, z1, nz) = 
    Grid(Float64, x0, x1, nx, y0, y1, ny, z0, z1, nz)

"""
    Grid([T=Float64,] x, y, z) -> ::Grid{T}

Create a new empty `Grid` by giving the coordinate ranges `x`, `y`
and `z`, optionally specifying the grid storage type `T`.
"""
function Grid(T, x::AbstractVector, y::AbstractVector, z::AbstractVector)
    x = map(T, x)
    y = map(T, y)
    z = map(T, z)
    Grid(x, y, z, Array{T}(undef, length(x), length(y), length(z)))
end
Grid(x, y, z) = Grid(Float64, x, y, z)

"""
    Grid(grid, data) -> ::Grid

Create a new `Grid` with the same dimensions and coordinates as `grid`,
but with different underlying `data`.
"""
function Grid(grid::Grid, data)
    size(grid) == size(data) ||
        throw(DimensionMismatch("grid and data array do not have the same dimensions"))
    Grid(grid.x, grid.y, grid.z, data)
end

Base.size(grid::Grid) = size(grid.data)
Base.getindex(grid::Grid, i::Int) = grid.data[i]
Base.getindex(grid::Grid,I::Vararg{Int,N}) where N = grid.data[I...]
Base.setindex!(grid::Grid, v, i::Int) = grid.data[i] = v
Base.setindex!(grid::Grid, v, I::Vararg{Int,N}) where N = grid.data[I...] = v
Base.axes(grid::Grid) = axes(grid.data)
Base.copy(grid::Grid) = Grid(grid.x, grid.y, grid.z, copy(grid.data))

function Base.show(io::IO, ::MIME"text/plain", grid::Grid{T}) where T
    print(io, """$(join(size(grid), '×')) Grid{$T}:
         x: $(grid.x)
         y: $(grid.y)
         z: $(grid.z)
         data: [values in range $(extrema(grid))]
        """)
end

"""
    MultiGrid{T,R} <: AbstractArray{T,3}

Set of several `Grid`s.  The first three dimensions are x, y and z,
and the fourth is sample.  E.g., to get the second sample in the chain,
do `grid[:,:,:,2]`.

# User-accessible fields
- `x`, `y` and `z`: Ranges containing the x, y and z coordinates in
  each dimension.  Units km.
- `grid`: 4D array containing the values at each point for each sample.
  Units km/s.

# Indexing
`MultiGrid`s are `AbstractArray`s, and so can be indexed just like normal
4D arrays with `getindex` and `setindex!` (e.g., `grid[:,:,3]`).
"""
struct MultiGrid{T,R} <: AbstractArray{T,4}
    x::R
    y::R
    z::R
    data::Array{T,4}
end

Base.size(grid::MultiGrid) = size(grid.data)
Base.getindex(grid::MultiGrid, i::Int) = grid.data[i]
Base.getindex(grid::MultiGrid ,I::Vararg{Int,N}) where N = grid.data[I...]
Base.setindex!(grid::MultiGrid, v, i::Int) = grid.data[i] = v
Base.setindex!(grid::MultiGrid, v, I::Vararg{Int,N}) where N = grid.data[I...] = v
Base.axes(grid::MultiGrid) = axes(grid.data)
Base.copy(grid::MultiGrid) = MultiGrid(grid.x, grid.y, grid.z, copy(grid.data))

"""
    MultiGrid(x0, x1, y0, y1, z0, z1, grid) -> ::MultiGrid

Construct a new `Grid` giving the start and end coordinates (`x0`
and `x1`) in each of the three dimensions, and a 3D array of values
`grid`.  The number of points in `grid` determines how many
coordinate points there are between e.g. `x0` and `x1`, and how
many differnt grids there are.
"""
function MultiGrid(x0, x1, y0, y1, z0, z1, data::AbstractArray{T,4}) where T
    nx, ny, nz, nt = size(grid)
    # Convert coordinates to correct precision
    x = map(T, range(x0, stop=x1, length=nx))
    y = map(T, range(y0, stop=y1, length=ny))
    z = map(T, range(z0, stop=z1, length=nz))
    MultiGrid(x, y, z, data)
end

"""
    MultiGrid(T, x0, x1, nx, y0, y1, ny, z0, z1, nz, nt) -> ::MultiGrid{T}
    MultiGrid(x0, x1, nx, y0, y1, ny, z0, z1, nz, nt) -> ::MultiGrid{Float64}

Create a new empty `Grid` by specifying the start and end values
of the coordinates in each dimension, the number of grid points, and
the number of grids (`nt`).
Specify the element type `T`, which defaults to `Float64` if not given.
"""
function MultiGrid(T, x0, x1, nx, y0, y1, ny, z0, z1, nz, nt)
    x = map(T, range(x0, stop=x1, length=nx))
    y = map(T, range(y0, stop=y1, length=ny))
    z = map(T, range(z0, stop=z1, length=nz))
    data = Array{T}(undef, nx, ny, nz, nt)    
    MultiGrid(x, y, z, data)
end
MultiGrid(x0, x1, nx, y0, y1, ny, z0, z1, nz, nt) = 
    MultiGrid(Float64, x0, x1, nx, y0, y1, ny, z0, z1, nz, nt)

"""
    MultiGrid(grid, data) -> ::Grid

Create a new `Grid` with the same dimensions and coordinates as `grid`,
but with different underlying `data`.
"""
function MultiGrid(grid::Grid, data)
    size(grid) == size(data) ||
        throw(DimensionMismatch("grid and data array do not have the same dimensions"))
    MultiGrid(grid.x, grid.y, grid.z, data)
end

"""
    MultiGrid(grids::AbstractArray{Grid})

Create a new `MultiGrid` from several `Grid`s.  All grids must
have the same coordinates, and therefore dimensions, and all
must have the same element type.
"""
function MultiGrid(grids::AbstractArray{<:Grid{T}}) where T
    isempty(grids) &&
        throw(ArgumentError("cannot construct a `MultiGrid`` from an empty set of `Grid`s"))
    all(==(coordinates(first(grids))), (coordinates(grid) for grid in grids)) ||
        throw(ArgumentError("not all grids have the same coordinates"))
    grid1 = first(grids)
    x, y, z = coordinates(grid1)
    nx, ny, nz = size(grid1)
    nt = length(grids)
    data = Array{T}(undef, nx, ny, nz, nt)
    for i in eachindex(grids)
        data[:,:,:,i] .= grids[i].data
    end
    MultiGrid(x, y, z, data)
end

#=
    Methods
=#

const GridOrMultiGrid = Union{Grid, MultiGrid}

"""
    coordinates(grid) -> (x, y, z)

Return a named tuple of the coordinates in each dimension of `grid`,
which may be a `Grid` or `MultiGrid`.
"""
coordinates(grid::GridOrMultiGrid) = (x=grid.x, y=grid.y, z=grid.z)

"""
    map_depth(f, grid) -> ::Vector

Perform the function `f` over all depth slices of the `grid` and
return the result.

# Example
```
julia> map(mean, grid)
```
"""
function map_depth(f::F, grid::GridOrMultiGrid) where F
    Base.map(f, eachslice(grid.data, dims=3))
end

"""
    map(f, grid::Union{Grid,MultiGrid}) -> grid′

Apply a function `f` to the underlying data of `grid`, and
return a new `Grid` with altered data.
"""
function Base.map(f::F, grid::GridOrMultiGrid) where F
    data = map(f, grid.data)
    MultiGrid(grid, data)
end

"""
    downsample(grid::GridOrMultiGrid, n) -> grid'

Create a new grid from `grid`, but taking every
`n`th point.  Note that this may mean that the
end coordinate of any dimension may not be preserved.
"""
function downsample(grid::Grid, n::Integer)
    x, y, z = coordinates(grid)
    Grid(x[begin:n:end], y[begin:n:end], z[begin:n:end],
        grid.data[begin:n:end,begin:n:end,begin:n:end])
end
function downsample(grid::MultiGrid, n::Integer)
    x, y, z = coordinates(grid)
    MultiGrid(x[begin:n:end], y[begin:n:end], z[begin:n:end],
        grid.data[begin:n:end,begin:n:end,begin:n:end,:])
end

end # module
