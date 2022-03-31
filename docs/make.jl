using QuantumESPRESSOFormatter
using Documenter

DocMeta.setdocmeta!(QuantumESPRESSOFormatter, :DocTestSetup, :(using QuantumESPRESSOFormatter); recursive=true)

makedocs(;
    modules=[QuantumESPRESSOFormatter],
    authors="Reno <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/QuantumESPRESSOFormatter.jl/blob/{commit}{path}#{line}",
    sitename="QuantumESPRESSOFormatter.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/QuantumESPRESSOFormatter.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/QuantumESPRESSOFormatter.jl",
)
