using MinkowskiReduction
using Documenter

DocMeta.setdocmeta!(MinkowskiReduction, :DocTestSetup, :(using MinkowskiReduction); recursive=true)

makedocs(;
    modules=[MinkowskiReduction],
    authors="Gus Hart",
    repo="https://github.com/glwhart/MinkowskiReduction.jl/blob/{commit}{path}#{line}",
    sitename="MinkowskiReduction.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://glwhart.github.io/MinkowskiReduction.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/glwhart/MinkowskiReduction.jl",
)
