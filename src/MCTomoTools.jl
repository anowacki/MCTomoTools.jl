"""
# MCTomoTools

Package for dealing with MCTomo runs.
"""
module MCTomoTools

export
    # Core types
    Chain,
    Model,
    
    # Helper module types
    Grid,
    MCTomoSettings,
    MultiGrid,

    # Core functions
    acceptance_ratios,
    sample_average_grid!,
    sample_average_grid,
    sample_chain!,
    sample_chain,
    sample_grid!,
    sample_grid,
    sample_source!,
    sample_source,
    sample_sources!,
    sample_sources,

    # Plotting
    plot_section,

    # Helper module functions
    # Accessors
    initial_model,
    # IO
    read_chains_burnin,
    read_initial_sample,
    read_receivers,
    read_raw_samples,
    read_settings,
    read_sources,
    read_temperatures,
    read_velocities!,
    read_velocities


using StaticArrays: FieldVector, SVector, @SVector

#=
    Modules not depending on anything else
=#
include("Grids.jl")
import .Grids

include("NameLists.jl")
import .NameLists

include("Settings.jl")
import .Settings

#=
    Core MCTomoTools functionality
=#
using .Grids: Grid
include("types.jl")
include("update.jl")
include("sample.jl")
include("qc.jl")

#=
    Helper modules
=#
include("Neighbours.jl")
import .Neighbours

include("InputOutput.jl")
import .InputOutput

include("Plot.jl")
import .Plot

# For reexport
using .InputOutput:
    read_chains_burnin,
    read_initial_sample,
    read_receivers,
    read_raw_samples,
    read_settings,
    read_sources,
    read_temperatures,
    read_velocities!,
    read_velocities
using .Settings: MCTomoSettings
using .Grids: Grid, MultiGrid
using .Plot: plot_section

end
