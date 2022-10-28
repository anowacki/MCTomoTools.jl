"""
# `Neighbours`

Module for dealing with nearest-neighbours, for the purposes of
computing the discretised Voronoi diagram from a set of vertices.
"""
module Neighbours

export
    fill_from_nodes!

using NearestNeighbors: KDTree, nn
using StaticArrays: SVector

using ..MCTomoTools: Nodes
import ..Grids
using ..Grids: Grid

"""
    fill_from_nodes!(grid::Grids.Grid, nodes::Nodes, property_number; method=:brute) -> grid

Fill a `grid` with the values of the `property_number`th property of a
set of `nodes`, constructing the discretised Voronoi tesselation of the
nodes.  In other words, for each point in `grid`, find the nearest node
in `nodes`, and fill that point with the value of `node.properties` in
position `property_number`.

By default (when `method` is `:brute`), a brute-force nearest-neighbour
search is performed, which does not allocate.  Pass `method=:kdtree` to
build a KD tree and query that in one go.
"""
function fill_from_nodes!(grid::Grid, nodes::Nodes{3,Nprop}, property_number;
    method=:brute
) where Nprop
    1 <= property_number <= Nprop ||
        throw(ArgumentError("requested property number is outside range"))
    
    if method === :brute
        # Using a brute-force method may be faster than creating the
        # KDTree and using that each time
        T = eltype(grid.x)
        
        for iz in eachindex(grid.z)
            z = grid.z[iz]
            for iy in eachindex(grid.y)
                y = grid.y[iy]
                for ix in eachindex(grid.x)
                    x = grid.x[ix]

                    mindist² = typemax(T)
                    nearest_inode = 0
                    for inode in eachindex(nodes.coords)
                        nx, ny, nz = T.(nodes.coords[inode])
                        dist² = (x - nx)^2 + (y - ny)^2 + (z - nz)^2
                        is_nearest = dist² < mindist²
                        mindist² = ifelse(is_nearest, dist², mindist²)
                        nearest_inode = ifelse(is_nearest, inode, nearest_inode)
                    end

                    grid[ix,iy,iz] = nodes.properties[nearest_inode][property_number]
                end
            end
        end
    elseif method === :kdtree
        # Method using the KDTree, created each time.  This is always
        # slower for the number of nodes we typically use (< 100).
        inds = _nearest_node_indices(grid, nodes)
        for i in eachindex(inds)
            @inbounds grid[i] = nodes.properties[inds[i]][property_number]
        end
    else
        throw(ArgumentError("method must be either `:brute` or `:kdtree`"))
    end

    grid
end

"""
    fill_from_nodes!(grids::Array{<:Grids.Grid}, nodes::Nodes) -> grids

As above, but fill in each grid within `grids`, where the first grid
gets the first property of each nodes, the second grid gets the second
property, and so on.
"""
function fill_from_nodes!(grids::AbstractArray{<:Grid}, nodes::Nodes{3,Nprop}) where Nprop
    length(grids) >= Nprop ||
        throw(ArgumentError("number of grids must be at least the number of properties"))
    all(g -> Grids.coordinates(g) == Grids.coordinates(first(grids)), grids) ||
        throw(ArgumentError("all grids must have the same coordinates"))
    inds = _nearest_node_indices(first(grids), nodes)
    for (igrid, iprop) in zip(eachindex(grids), 1:Nprop)
        grid = grids[igrid]
        for (i, ind) in enumerate(inds)
            grid[i] = nodes.properties[ind][iprop]
        end
    end
    grids
end

"""
    _nearest_node_indices(grid, nodes) -> indices

For each point in `grid`, calculate the index of the nearest node in 
`nodes` and return this as a vector.  Therefore, the closest point to
`grid[i]` is `nodes[indices[i]]`.
"""
function _nearest_node_indices(grid::Grid, nodes)
    # The nodes may be different each time, so create a new tree each time
    tree = KDTree(nodes.coords)

    # Do the lookup on all grid points at once as this is usually quicker
    grid_coords = Iterators.product(Grids.coordinates(grid)...)
    inds, distances = nn(tree, vec([SVector(coord) for coord in grid_coords]))
    inds
end

"""
_nearest_node_indices(coords, nodes) -> indices

Return the `indices` of the nearest node to each coordinate in
`coords`.  This should be a `Vector{SVector{3}}` of the coordinates
of each point of interest.
"""
function _nearest_node_indices(coords, nodes)
    # Create a new tree each time as the nodes may have moved
    tree = KDTree(nodes.coords)
    inds, distances = nn(tree, coords)
    inds
end

"""
    _grid_coordinates(grid) -> coordinates::Vector{SVector{3}}

Calculate a vector of coordinates for every point in `grid`.
This may be passed to `_nearest_node_indices` rather than the
`Grid` itself to avoid repeatedly allocating this vector of
coordinates.
"""
function _grid_coordinates(grid::Grid)
    coords = Iterators.product(Grids.coordinates(grid)...)
    vec([SVector(coord) for coord in grid_coords])
end

end # module
