module PWscf

using Compat: eachrow
using Crystallography: ReciprocalPoint, MonkhorstPackGrid
using Formatting: sprintf1
using QuantumESPRESSOBase.Inputs.PWscf:
    AtomicSpecies,
    AtomicPosition,
    ReciprocalPoint,
    AtomicSpeciesCard,
    AtomicPositionsCard,
    CellParametersCard,
    AtomicForce,
    SpecialPointsCard,
    GammaPointCard,
    KMeshCard,
    PWInput,
    optionof

import AbInitioSoftwareBase.Inputs: FormatConfig, asstring

export format_file, format_text

FormatConfig(
    ::Union{
        AtomicSpecies,
        AtomicPosition,
        ReciprocalPoint,
        AtomicSpeciesCard,
        AtomicPositionsCard,
        CellParametersCard,
        MonkhorstPackGrid,
        AtomicForce,
        SpecialPointsCard,
        GammaPointCard,
        KMeshCard,
    },
) = FormatConfig(;
    delimiter = " ",
    newline = "\n",
    indent = ' '^4,
    float = "%14.9f",
    int = "%5d",
    bool = ".%."
)

"""
    asstring(data::AtomicSpecies)

Return a `String` representing a `AtomicSpecies`, valid for Quantum ESPRESSO's input.
"""
function Base.print(io::IO, data::AtomicSpecies)
    config = FormatConfig(data)
    print(io,
        join(
            (
                config.indent,
                sprintf1("%3s", data.atom),
                sprintf1(config.float, data.mass),
                data.pseudopot,
            ),
            config.delimiter,
        )
    )
    return nothing
end
"""
    asstring(card::AtomicSpeciesCard)

Return a `String` representing a `AtomicSpeciesCard`, valid for Quantum ESPRESSO's input.
"""
function Base.print(io::IO, card::AtomicSpeciesCard)
    config = FormatConfig(card)
    data = union(card.data)
    print(io, join(("ATOMIC_SPECIES", map(string, data)...), config.newline))
    return nothing
end
"""
    asstring(data::AtomicPosition)

Return a `String` representing a `AtomicPosition`, valid for Quantum ESPRESSO's input.
"""
function Base.print(io::IO, data::AtomicPosition)
    config = FormatConfig(data)
    content = join(
        (
            config.indent,
            sprintf1("%3s", data.atom),
            map(x -> sprintf1(config.float, x), data.pos)...,
        ),
        config.delimiter,
    )
    if !all(data.if_pos)
        print(io, join((content, map(Int, data.if_pos)...), config.delimiter))
    else
        print(io, content)
    end
    return nothing
end
"""
    asstring(card::AtomicPositionsCard)

Return a `String` representing a `AtomicPositionsCard`, valid for Quantum ESPRESSO's input.
"""
function Base.print(io::IO, card::AtomicPositionsCard)
    config = FormatConfig(card)
    print(io, join(
        ("ATOMIC_POSITIONS { $(optionof(card)) }", map(string, card.data)...),
        config.newline,
    ))
    return nothing
end
"""
    asstring(card::CellParametersCard)

Return a `String` representing a `CellParametersCard`, valid for Quantum ESPRESSO's input.
"""
function Base.print(io::IO, card::CellParametersCard)
    config = FormatConfig(card)
    print(io, join(
        (
            "CELL_PARAMETERS { $(optionof(card)) }",
            map(eachrow(card.data)) do row
                join((sprintf1(config.float, x) for x in row))
            end...,
        ),
        config.newline,
    ))
    return nothing
end
"""
    asstring(data::MonkhorstPackGrid)

Return a `String` representing a `MonkhorstPackGrid`, valid for Quantum ESPRESSO's input.
"""
function Base.print(io::IO, data::MonkhorstPackGrid)
    config = FormatConfig(data)
    print(io, config.indent * join(map([data.mesh; data.is_shift]) do x
            sprintf1(config.int, x)
        end, config.delimiter))
    return nothing
end
"""
    asstring(data::SpecialKPoint)

Return a `String` representing a `SpecialKPoint`, valid for Quantum ESPRESSO's input.
"""
function Base.print(io::IO, data::ReciprocalPoint)
    config = FormatConfig(data)
    print(io, config.indent * join(
        map(x -> sprintf1(config.float, x), [data.coord..., data.weight]),
        config.delimiter,
    ))
    return nothing
end
"""
    asstring(card::KPointsCard)

Return a `String` representing a `KPointsCard`, valid for Quantum ESPRESSO's input.
"""
function Base.print(io::IO, card::SpecialPointsCard)
    config = FormatConfig(card)
    content = "K_POINTS { $(optionof(card)) }" * config.newline
    print(io, join((content, length(card.data), map(string, card.data)...), config.newline))
    return nothing
end
function Base.print(io::IO, card::GammaPointCard)
    config = FormatConfig(card)
    print(io, "K_POINTS { $(optionof(card)) }" * config.newline)
    return nothing
end
function Base.print(io::IO, card::KMeshCard)
    config = FormatConfig(card)
    content = "K_POINTS { $(optionof(card)) }" * config.newline
    print(io, content * string(card.data))
    return nothing
end

function format_file(filename::AbstractString; overwrite::Bool = true, kwargs...)
    text = read(filename, String)
    formatted_text = format_text(text; kwargs...)
    if overwrite
        open(filename, "w") do io
            write(io, formatted_text)
        end
    else
        println(formatted_text)
    end
end

format_text(text::AbstractString) = string(parse(PWInput, text))

end
