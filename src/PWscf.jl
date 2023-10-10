module PWscf

using CrystallographyBase: MonkhorstPackGrid
using Formatting: sprintf1
using QuantumESPRESSOBase.PWscf:
    SpecialPoint,
    AtomicSpecies,
    AtomicPosition,
    AtomicSpeciesCard,
    AtomicPositionsCard,
    CellParametersCard,
    AtomicForce,
    SpecialPointsCard,
    GammaPointCard,
    KMeshCard,
    PWInput,
    getoption
using QuantumESPRESSOParser.PWscf

import AbInitioSoftwareBase: FormatConfig

export format_file, format_text

FormatConfig(
    ::Union{
        AtomicSpecies,
        AtomicPosition,
        SpecialPoint,
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
    delimiter=" ", newline="\n", indent=' '^4, float="%14.9f", int="%5d", bool=".%."
)

function Base.print(io::IO, data::AtomicSpecies)
    config = FormatConfig(data)
    print(
        io,
        join(
            (
                config.indent,
                sprintf1("%3s", data.atom),
                sprintf1(config.float, data.mass),
                data.pseudopot,
            ),
            config.delimiter,
        ),
    )
    return nothing
end
function Base.print(io::IO, card::AtomicSpeciesCard)
    config = FormatConfig(card)
    data = union(card.data)
    print(io, join(("ATOMIC_SPECIES", map(string, data)...), config.newline))
    return nothing
end
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
function Base.print(io::IO, card::AtomicPositionsCard)
    config = FormatConfig(card)
    print(
        io,
        join(
            ("ATOMIC_POSITIONS { $(getoption(card)) }", map(string, card.data)...),
            config.newline,
        ),
    )
    return nothing
end
function Base.print(io::IO, card::CellParametersCard)
    config = FormatConfig(card)
    print(
        io,
        join(
            (
                "CELL_PARAMETERS { $(getoption(card)) }",
                Iterators.map(eachrow(card.data)) do row
                    join((sprintf1(config.float, x) for x in row))
                end...,
            ),
            config.newline,
        ),
    )
    return nothing
end
function Base.print(io::IO, data::MonkhorstPackGrid)
    config = FormatConfig(data)
    print(
        io,
        config.indent * join(
            Iterators.map([data.mesh; data.is_shift]) do x
                sprintf1(config.int, x)
            end,
            config.delimiter,
        ),
    )
    return nothing
end
function Base.print(io::IO, data::SpecialPoint)
    config = FormatConfig(data)
    print(
        io,
        config.indent * join(
            Iterators.map(x -> sprintf1(config.float, x), [data.coord..., data.weight]),
            config.delimiter,
        ),
    )
    return nothing
end
function Base.print(io::IO, card::SpecialPointsCard)
    config = FormatConfig(card)
    content = "K_POINTS { $(getoption(card)) }" * config.newline
    print(io, join((content, length(card.data), map(string, card.data)...), config.newline))
    return nothing
end
function Base.print(io::IO, card::GammaPointCard)
    config = FormatConfig(card)
    print(io, "K_POINTS { $(getoption(card)) }" * config.newline)
    return nothing
end
function Base.print(io::IO, card::KMeshCard)
    config = FormatConfig(card)
    content = "K_POINTS { $(getoption(card)) }" * config.newline
    print(io, content * string(card.data))
    return nothing
end

function format_file(filename::AbstractString; overwrite::Bool=true, kwargs...)
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
