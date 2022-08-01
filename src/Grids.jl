"""
# `Grids`

Module for dealing with velocity grids
"""
module Grids

export
    Grid,
    coordinates

"Alias for the type created by `range` or `start:step:end`"
const _Range{T} = Base.StepRangeLen{T, Base.TwicePrecision{T}, Base.TwicePrecision{T}, Int}

"""
    Grid{T} <: AbstractArray{T,3}

Grid containing a 3D volume of points spaced equally.

# User-accessible fields
- `x`, `y` and `z`: Ranges containing the x, y and z coordinates in
  each dimension.  Units km.
- `grid`: 3D array containing the values at each point.  Units km/s.

# Indexing
`Grid`s are `AbstractArray`s, and so can be indexed just like normal
3D arrays with `getindex` and `setindex!` (e.g., `grid[:,:,3]`).
"""
struct Grid{T} <: AbstractArray{T,3}
    x::_Range{T}
    y::_Range{T}
    z::_Range{T}
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

"""
    coordinates(grid) -> (x, y, z)

Return a named tuple of the coordinates in each dimension of `grid`.
"""
coordinates(grid::Grid) = (x=grid.x, y=grid.y, z=grid.z)

"""
    map_depth(f, grid) -> ::Vector

Perform the function `f` over all depth slices of the `grid` and
return the result.

# Example
```
julia> map(mean, grid)
```
"""
function map_depth(f::F, grid::Grid) where F
    Base.map(f, eachslice(grid.data, dims=3))
end

"""
    map(f, grid::Grid) -> ::Grid

Apply a function `f` to the underlying data of `grid`, and
return a new `Grid` with altered data.
"""
function Base.map(f::F, grid::Grid) where F
    data = map(f, grid.data)
    Grid(grid, data)
end

end # module
