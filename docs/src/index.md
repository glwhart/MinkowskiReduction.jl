```@meta
CurrentModule = MinkowskiReduction
DocTestSetup = quote
    using MinkowskiReduction
    using LinearAlgebra
end
```

# MinkowskiReduction.jl

A Julia package for Minkowski reduction of three-dimensional lattice
bases. Given any basis for a 3D lattice, `minkReduce` returns the
basis with the shortest possible vectors — equivalently, the basis
closest to orthogonal.

Used in crystallography and density-functional-theory calculations to
normalise a crystal's unit cell before further analysis (e.g.
[Spacey.jl](https://github.com/glwhart/Spacey.jl) uses it as the
first step of spacegroup detection).

## Quick taste

```jldoctest
julia> using MinkowskiReduction

julia> M = [1.0  0.0  1.0;
            0.0  1.0  1.0;
            0.0  0.0  1.0];     # a skewed cubic basis

julia> R, P = minkReduce(M);

julia> R                         # reduced basis (orthonormal)
3×3 Matrix{Float64}:
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0

julia> M * P == R                # P is the exact integer transform
true
```

## The documentation is organised by what you need

This documentation follows the
[Diátaxis](https://diataxis.fr/) framework. Four kinds of pages serve
four different needs:

- **[Tutorial](tutorial.md)** — *learning*. A single guided walk-through
  that takes you from installation to a working reduction. Start here
  if you're new.

- **[How-to guides](how-to.md)** — *getting things done*. Short recipes
  for specific tasks: transform atomic positions, check reducedness,
  reduce a 2D sublattice, etc. Use these when you know what you want
  to do.

- **Reference** — *looking things up*.
  - [API reference](reference/api.md) — every exported name, with full
    signatures.
  - [The 12 Minkowski conditions](reference/conditions.md) — the
    mathematical contract that defines reducedness.

- **Explanation** — *understanding*.
  - [The algorithm](explanation/algorithm.md) — how `minkReduce`
    actually works (the Nguyen–Stehlé greedy algorithm).
  - [Non-uniqueness](explanation/non-uniqueness.md) — why the same
    lattice can produce different reduced bases.
  - [Precision](explanation/precision.md) — floating-point behaviour,
    scale-aware tolerances, and the `Int64`-overflow pitfall.
  - [Context](explanation/context.md) — how Minkowski reduction
    compares to LLL, Niggli, and Selling reduction, and why this
    package exists.

## Installation

```julia
julia> using Pkg

julia> Pkg.add("MinkowskiReduction")
```

The package has no external dependencies — only `LinearAlgebra`,
`Random`, and `Test` from the standard library.
