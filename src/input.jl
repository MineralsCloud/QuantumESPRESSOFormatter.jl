using IterTools: imap
using PyFortran90Namelists: fstring
using QuantumESPRESSOBase: QuantumESPRESSOInput, dropdefault, groupname

import AbInitioSoftwareBase.Inputs: FormatConfig, Namelist

FormatConfig(::Union{QuantumESPRESSOInput,Namelist}) = FormatConfig(;
    delimiter=" ", newline="\n", indent=' '^4, float="%f", int="%i", bool=".%."
)

function Base.print(io::IO, input::QuantumESPRESSOInput)
    newline = FormatConfig(input).newline
    iter = imap(1:nfields(input)) do i
        x = getfield(input, i)
        x === nothing ? "" : string(x)
    end
    println(io, join(iter, newline) * newline)  # Add a new line at the end of line to prevent errors
    return nothing
end
function Base.print(io::IO, nml::Namelist)
    dict = dropdefault(nml)
    config = FormatConfig(nml)
    indent, delimiter, newline = config.indent, config.delimiter, config.newline
    iter = imap(dict) do (key, value)
        if value isa AbstractVector
            data = imap(enumerate(value)) do (i, x)
                if x !== nothing
                    indent * join((string(key, '(', i, ')'), "=", fstring(x)), delimiter)
                end
            end
            join(Iterators.filter(!isnothing, data), newline)
        elseif value isa NamedTuple
            data = imap(value) do (x, y)
                indent * join((string(key, '%', x), "=", fstring(y)), delimiter)
            end
            join(data, newline)
        else
            indent * join((string(key), "=", fstring(value)), delimiter)
        end
    end
    content = join(iter, newline)
    print(io, join(filter(!isempty, ("&" * groupname(nml), content, '/')), newline))
    return nothing
end
