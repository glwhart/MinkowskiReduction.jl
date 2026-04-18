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
    repo=Documenter.Remotes.GitHub("glwhart", "MinkowskiReduction.jl"),
    sitename="MinkowskiReduction.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://glwhart.github.io/MinkowskiReduction.jl",
        assets=String[],
    ),
    pages=[
        "Home"        => "index.md",
        "Tutorial"    => "tutorial.md",
        "How-to"      => "how-to.md",
        "Reference"   => [
            "reference/api.md",
            "reference/conditions.md",
        ],
        "Explanation" => [
            "explanation/algorithm.md",
            "explanation/non-uniqueness.md",
            "explanation/precision.md",
            "explanation/context.md",
        ],
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
