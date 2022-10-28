module Settings

export getpath

"""
    MCTomoSettings

Structure holding a set of settings for an MCTomo run, typically in
a file called 'MCTomo.inp'.

The structure also holds the directory
"""
struct MCTomoSettings
    settings::Dict{String,Dict{String,Dict{String,Any}}}
    directory::String
end

"""
    getindex(settings::MCTomoSettings, section, struct, field)

Get the value from an `MCTomoSettings` structure, where `section`
is the name of the namelist section, `struct` is the name of
the struct which is read in, and `field` is the name of the field
within that derived type.
"""
Base.getindex(settings::MCTomoSettings, section, struc, field) =
    settings.settings[section][struc][field]

"""
    joinpath(settings::MCTomoSettings, args...) -> path

Join `args` to the base directory of an MCTomo run, as held in a
`MCTomoSettings` object.
"""
Base.joinpath(settings::MCTomoSettings, args...) =
    joinpath(settings.directory, args...)

"""
    getpath(settings, section, struc, field) -> path

Return the relative path to a file stored in a `MCTomoSettings` file
`settings`.
"""
function getpath(settings::MCTomoSettings, section, struc, field)
    value = settings[section, struc, field]
    value isa String ||
        throw(ArgumentError("&$section $struc%$field is not a string value (is $value)"))
    joinpath(settings.directory, value)
end

end # module
