"""
# `NameLists`

Module for reading Fortran namelist files.

!!! note
    This module only supports reading derived types and their
    elements, not 'bare' variables, arrays, or many other things.
    It is only meant to be used with MCTomoTools!

## Functions
- [`read_namelist`](@ref)
"""
module NameLists

export read_namelist

"Order in which data are parsed, finally returning a string"
TYPE_ORDER = (Bool, Int, Float64, Complex{Float64}, String)

struct Line
    new_section::Union{String,Nothing}
    end_of_section::Bool
    type::Union{String,Nothing}
    field::Union{String,Nothing}
    value::Union{String,Nothing}
    empty::Bool
end
Line(; new_section=nothing, end_of_section=false, type=nothing, field=nothing, value=nothing, empty=false) =
    Line(new_section, end_of_section, type, field, value, empty)

"""
    read_namelist(file_or_io; types=Dict()) -> ::Dict{String,Dict{String,Any}}

Read a set of values from a Fortran name list file.  (See [`NameLists`](@ref)
for an example of the structue.)  All section, variable and field names,
and hence all keys of the `Dict`s, are converted to lowercase, as case
is not important in namelist files.

The current implementation only supports the reading in of values into
derived types!  Reading in of 'bare' variables is not supported.

This function will attempt to read data as types in the following order:

1. Logical (Julia type `Bool`)
2. Integer (`Int`)
3. Real (`Float64`)
4. Complex (`Complex{Float64}`)
5. Character (`String`) (but strings should be enclosed in quotes in
   any case and this should be automatically recognised)

The exception to this is if the user supplies a `Dict` which maps
the section, variable and field names onto a specific type, in which
case those values are attempted to be read as that type, and an error
is thrown if they cannot be.
"""
function read_namelist(io::IO;
    types=Dict{String,Dict{String,Dict{String,DataType}}}(),
    file=nothing
)
    namelist = Dict{String,Dict{String,Dict{String,Any}}}()

    local current_section

    for (linenumber, linestring) in enumerate(eachline(io))
        line = _get_line(linestring, linenumber, file)
        (line.empty || line.end_of_section) && continue

        if !isnothing(line.new_section)
            current_section = line.new_section
            namelist[current_section] = Dict{String,Dict{String,Any}}()
        else
            if !haskey(namelist[current_section], line.type)
                namelist[current_section][line.type] = Dict{String,Any}()
            end
            value = _parse_value(current_section, line, types, linenumber, file)
            namelist[current_section][line.type][line.field] = value
        end
    end
    namelist
end

read_namelist(file; kwargs...) =
    open(io -> read_namelist(io; file=file, kwargs...), file)

"""
Read the `line` and return a `Line` with fields `new_section::String`,
`end_of_section::Bool`, `type::String`, `field::String` and `value::String`,
`empty::Bool`, where all fields may also be `nothing`.

If the current line starts a new section, then only `new_section` is not
`nothing`; if this is a line with a type and field, then `type`, `field` and `value` are
not `nothing`.  An end of section line has `end_of_section = true`.  An empty
line (with optional comments) has `empty = true`.

!!! note
    The current implementation does not support strings which include
    the `!` character.
"""
function _get_line(line, linenumber=nothing, file=nothing)
    # Empty or all comment lines
    (isempty(line) || occursin(r"^\s*!", line) || occursin(r"^\s*$", line)) &&
        return Line(; empty=true)

    # Remove comments and leading/trailing whitespace
    line = strip(replace(line, r"!.*"=>""))

    # New section
    if occursin(r"^&", line)
        tokens = split(line, r"[&!]")
        length(tokens) >= 2 || error("too few tokens to parse a section name in '$line'" *
            _line_number_string(linenumber, file))
        new_section = lowercase(strip(tokens[2]))
        return Line(; new_section)

    # End of section
    elseif occursin(r"^/(([\s]*!+.*$)|([\s!]*$))", line)
        return Line(; end_of_section=true)

    # Should be a line with values on, maybe with a trailing comma
    else
        matches = match(r"([A-Za-z_][A-Za-z0-9_]*)%([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([A-Za-z0-9\\._+\-\\'\"]+)\s*,?$",
            line)
        (isnothing(matches) || length(matches) < 3) &&
            error("cannot get names and values from '$line'" *
            _line_number_string(linenumber, file))
        return Line(; type=lowercase(matches[1]), field=lowercase(matches[2]),
            value=matches[3])
    end
end

function _parse_value(current_section, line, types, linenumber=nothing, file=nothing)
    requested_type = if haskey(types, current_section)
        if haskey(types[current_section], line.type)
            if haskey(types[current_section][line.type], line.field)
                types[current_section][line.type][line.field]
            else
                nothing
            end
        else
            nothing
        end
    else
        nothing
    end

    if !isnothing(requested_type)
        value = _tryparse(requested_type, line.value)
        isnothing(value) &&
            error("cannot parse '$(line.value)' as $requested_type" *
                _line_number_string(linenumber, file))
        return value
    else
        for type in TYPE_ORDER
            value = _tryparse(type, line.value)
            if !isnothing(value)
                return value
            end
        end
        error("do not know how to parse '$(line.value)'" *
            _line_number_string(linenumber, file))
    end
end

function _tryparse(::Type{String}, token)
    if occursin(r"^(\'.*\'|\".*\")$", token)
        # This is okay because we have established the first and last
        # characters are either ' or "
        token[begin+1:end-1]
    else
        nothing
    end
end
function _tryparse(::Type{Bool}, token)
    if startswith(token, ".T") || token == "T"
        true
    elseif startswith(token, ".F") || token == "F"
        false
    else
        nothing
    end
end
_tryparse(T, token) = tryparse(T, token)

function _line_number_string(linenumber, file)
    if isnothing(linenumber)
        ""
    else
        if isnothing(file)
            " (line $linenumber)"
        else
            " (line $linenumber of file '$file')"
        end
    end
end

end # module
