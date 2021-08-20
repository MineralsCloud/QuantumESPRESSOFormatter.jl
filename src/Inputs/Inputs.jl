module Inputs

using Compat
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
    iter = Iterators.map(1:nfields(input)) do i
        x = getfield(input, i)
        x === nothing ? "" : asstring(x)
    end
    return join(iter, newline) * newline  # Add a new line at the end of line to prevent errors
end
"""
    asstring(nml::Namelist)

Return a `String` representing a `Namelist`, valid for Quantum ESPRESSO's input.
"""
function asstring(nml::Namelist)
    dict = dropdefault(nml)
    config = FormatConfig(nml)
    indent, delimiter, newline = config.indent, config.delimiter, config.newline
    iter = Iterators.map(dict) do (key, value)
        if value isa AbstractVector
            data = Iterators.map(enumerate(value)) do (i, x)
                if x !== nothing
                    indent * join((string(key, '(', i, ')'), "=", fstring(x)), delimiter)
                end
            end
            join(data, newline)
        elseif value isa NamedTuple
            data = Iterators.map(value) do (x, y)
                indent * join((string(key, '%', x), "=", fstring(y)), delimiter)
            end
            join(data, newline)
        else
            indent * join((string(key), "=", fstring(value)), delimiter)
        end
    end
    content = join(iter, newline)
    return join(filter(!isempty, ("&" * groupname(nml), content, '/')), newline)
end

include("PWscf.jl")

end
