using QuantumESPRESSOFormatter
using Documenter

DocMeta.setdocmeta!(QuantumESPRESSOFormatter, :DocTestSetup, :(using QuantumESPRESSOFormatter); recursive=true)

makedocs(;
    modules=[QuantumESPRESSOFormatter],
    authors="singularitti <singularitti@outlook.com> and contributors",
    repo="https://github.com/MineralsCloud/QuantumESPRESSOFormatter.jl/blob/{commit}{path}#{line}",
    sitename="QuantumESPRESSOFormatter.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/QuantumESPRESSOFormatter.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Installation Guide" => "installation.md",
        ],
        # "Public API" => "public.md",
        "Developer Docs" => [
            "Contributing" => "developers/contributing.md",
            "Style Guide" => "developers/style-guide.md",
            "Design Principles" => "developers/design-principles.md",
        ],
        "Troubleshooting" => "troubleshooting.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/QuantumESPRESSOFormatter.jl",
    devbranch="main",
)
