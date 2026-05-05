```@meta
CurrentModule = MinkowskiReduction
DocTestSetup = quote
    using MinkowskiReduction
    using LinearAlgebra
end
```

# MinkowskiReduction.jl

A Julia package for Minkowski reduction of three-dimensional lattice
bases (or Gauss reduction of a two-dimensional lattice). Given any basis for a 3D lattice, `mink_reduce` returns the
basis with the shortest possible vectors — equivalently, the basis
closest to orthogonal.

Minkowski reduction is useful in crystallography and density-functional-theory calculations to put a crystal's unit cell into a canonical form before further analysis (e.g.
[Spacey.jl](https://github.com/glwhart/Spacey.jl) uses it as the
first step of spacegroup detection).

## Quick taste

```jldoctest
julia> using MinkowskiReduction

julia> M = [1.0  0.0  1.0;
            0.0  1.0  1.0;
            0.0  0.0  1.0];     # a skewed cubic basis

julia> R, P = mink_reduce(M);

julia> R                         # reduced basis (orthonormal)
3×3 Matrix{Float64}:
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0

julia> M * P == R                # P is the exact integer transform taking the original basis to the reduced one
true
```

## The documentation is organised by what you need

This documentation follows the
[Diátaxis](https://diataxis.fr/) framework. Four kinds of pages serve
four different kinds of user:

- **[Tutorial](tutorial.md)** — *learning*. A guided walk-through
  that takes you from installation to reducing a lattice. Start here
  if you're new.

- **[How-to guides](how-to.md)** — *getting things done*. Short recipes
  for specific tasks: reduce a lattice basis (3D `mink_reduce` and 2D `gauss_reduce`), check the orthogonality defect of a lattice `orthogonality_defect`, check if a lattice basis is already reduced `is_mink_reduced`, etc. Use these when you know what you want
  to do.

- **Reference** — *looking things up*. Use this if you just need to know how to call a function.
  - [API reference](reference/api.md) — every exported name, with full
    signatures.
  - [The 12 Minkowski conditions](reference/conditions.md) — the
    mathematical contract that defines Minkowski reduction (shortest basis vectors, most orthogonal basis).

- **Explanation** — *understanding*. The math and logic behind the algorithm.
  - [The algorithm](explanation/algorithm.md) — how `mink_reduce`
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
