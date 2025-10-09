# Plotting function stubs, extended by MCTomoToolsMakieExt when
# Makie is loaded

"""
    plot_section_makie([gridposition,] chain::Union{Chain,AbstractArray{<:Chain}}, grid::Grid, property::Symbol; source=nothing, receivers=nothing, burnin=0, thin=1000, x=nothing, y=nothing, z=nothing, method=:brute) -> ::Makie.Figure

Create a cross-section showing the mean value of `property` for every
`thin`th sample, starting at sample `burnin + 1`, from the `chain`
of samples.  Specify where the cross-section is at either, `x`, `y` or
`z` km.

# Keyword arguments
- `sources`: Named tuple or struct where fields `x`, `y` and `z` respectively
  give the x, y and z coordinates (km) of the body-wave sources
- `receivers`: Named tuple or struct where fields `x`, `y` and `z` respectively
  give the x, y and z coordinates (km) of the receivers
"""
function plot_section_makie end
