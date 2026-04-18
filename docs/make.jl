using MinkowskiReduction
using Documenter
using LinearAlgebra  # `det` is referenced by several docstring examples

DocMeta.setdocmeta!(
    MinkowskiReduction,
    :DocTestSetup,
    :(using MinkowskiReduction; using LinearAlgebra);
    recursive=true,
)

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
    doctest=true,   # explicit: fail the build if any jldoctest example is stale
)

deploydocs(;
    repo="github.com/glwhart/MinkowskiReduction.jl", 
    devbranch = "main",
    devurl="dev",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#" ]
)
