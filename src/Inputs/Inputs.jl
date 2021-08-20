module Inputs

using PyFortran90Namelists: fstring
using QuantumESPRESSOBase.Inputs: QuantumESPRESSOInput, dropdefault, groupname

import AbInitioSoftwareBase.Inputs: FormatConfig, Namelist, asstring

FormatConfig(::Union{QuantumESPRESSOInput,Namelist}) = FormatConfig(;
    delimiter = " ",
    newline = "\n",
    indent = ' '^4,
    float = "%f",
    int = "%i",
    bool = ".%.",
)

"""
    asstring(input::QuantumESPRESSOInput)

Return a `String` representing a `QuantumESPRESSOInput`, valid for Quantum ESPRESSO's input.
"""
function asstring(input::QuantumESPRESSOInput)
    newline = FormatConfig(input).newline
    return join(map(_asstring, getfield(input, i) for i in 1:nfields(input)), newline) *
           newline  # Add a new line at the end of line to prevent errors
end
"""
    asstring(nml::Namelist)

Return a `String` representing a `Namelist`, valid for Quantum ESPRESSO's input.
"""
function asstring(nml::Namelist)
    dict = dropdefault(nml)
    config = FormatConfig(nml)
    indent, delimiter, newline = config.indent, config.delimiter, config.newline
    iter = (
        if value isa AbstractVector
            data = (
                indent * join((string(key, '(', i, ')'), "=", fstring(x)), delimiter) for
                (i, x) in enumerate(value) if !isnothing(x)
            )
            join(data, newline)
        elseif value isa NamedTuple
            data = (
                indent * join((string(key, '%', x), "=", fstring(y)), delimiter) for
                (x, y) in value
            )
            join(data, newline)
        else
            indent * join((string(key), "=", fstring(value)), delimiter)
        end for (key, value) in dict
    )
    content = join(iter, newline)
    return join(filter(!isempty, ("&" * groupname(nml), content, '/')), newline)
end
_asstring(::Nothing) = ""
_asstring(x) = asstring(x)

include("PWscf.jl")

end
