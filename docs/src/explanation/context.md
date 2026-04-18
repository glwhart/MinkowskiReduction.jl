```@meta
CurrentModule = MinkowskiReduction
```

# Why Minkowski reduction?

A single lattice `L ⊂ ℝ³` has *infinitely many* bases — any two bases
related by a unimodular integer matrix span the same lattice. Some of
these bases are pleasant to work with (short vectors, near-orthogonal);
most are not. A **reduction** is a procedure that takes an arbitrary
basis for `L` and replaces it with one of the pleasant ones.

Several reduction notions compete. This page explains where Minkowski
reduction sits among them, and why this package chose it.

## What Minkowski reduction is

A basis `{b₁, b₂, b₃}` is **Minkowski reduced** if, for each `i`,
`bᵢ` is the shortest lattice vector such that `{b₁, …, bᵢ}` can be
extended to a full basis of `L`. Equivalently — and this is the
result that makes low-dimensional Minkowski reduction attractive —
the norms of a Minkowski-reduced basis are exactly the **successive
minima** `λ₁, λ₂, λ₃` of the lattice:

- `λ₁` is the length of the shortest non-zero lattice vector,
- `λ₂` is the length of the shortest vector linearly independent
  from the shortest,
- `λ₃` is the length of the shortest vector linearly independent
  from the first two.

So a Minkowski-reduced basis achieves *all three* successive minima
simultaneously. It is, in a precise sense, the shortest possible
basis for the lattice.

This theorem holds in dimensions ≤ 4. In dimensions ≥ 5 the
successive minima may not be attainable by any single basis, and the
greedy algorithm of [Nguyen & Stehlé](algorithm.md) can fail to even
find `λ₁`. Fortunately for crystallographers, physical lattices are
three-dimensional.

## Orthogonality defect

A related invariant of a basis is its **orthogonality defect**:

```
δ(b₁, b₂, b₃) = ‖b₁‖ · ‖b₂‖ · ‖b₃‖ / |det|
```

By Hadamard's inequality, `δ ≥ 1`, with equality exactly when the
basis is orthogonal. For skewed bases, `δ` grows without bound as
the basis becomes more ill-conditioned. A Minkowski-reduced basis
minimises `δ` over all bases of the lattice.

[`orthogonalityDefect`](@ref) computes this. Representative values:

- Simple cubic primitive cell: `δ = 1` exactly.
- FCC primitive (`a₁₂₃ = {(0,½,½), (½,0,½), (½,½,0)}`): `δ = √2 ≈ 1.414`.
- BCC primitive: `δ = 3√3/4 ≈ 1.299`.
- Hexagonal primitive: `δ = 2/√3 ≈ 1.155`.

These are what a correct reducer **must** achieve starting from any
basis for the corresponding lattice. The test suite verifies this for
all seven of these lattice types plus a rhombohedral example.

## Comparison to related reductions

### LLL reduction

[LLL (Lenstra–Lenstra–Lovász, 1982)](https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm)
is the famous polynomial-time reduction that works in arbitrary
dimension. Its first basis vector satisfies

```
‖b₁‖ ≤ 2^{(d-1)/2} · λ₁
```

— an **exponential approximation factor** in the dimension `d`.
This is vastly weaker than Minkowski reduction's exact achievement
of `λ₁`, but LLL's tractability across all `d` makes it the tool of
choice for cryptography and integer programming. In 3D, Minkowski
reduction (via the greedy algorithm) is strictly better: it achieves
exact optimality *and* has quadratic bit-complexity.

### Niggli reduction

[Niggli reduction](https://dictionary.iucr.org/Niggli_reduced_cell) is
a 3D crystallographic standard that adds sign and axis normalisations
on top of Buerger reduction (itself equivalent to Minkowski reduction
in 3D) so that the final basis is a **unique function of the lattice**.
This is what crystallography databases use to canonicalise a cell for
lookup.

This package deliberately stops at Minkowski reduction — see
[Non-uniqueness](non-uniqueness.md). If you need a canonical form for
cell comparison, post-process with a Niggli reducer. If you only need
point-group or spacegroup analysis (which only depends on the
lattice's abstract symmetry), the non-canonicalness doesn't matter.

### Selling reduction

[Selling reduction](https://doi.org/10.1107/S0108767395001268) is a
different 3D approach that works on four vectors
`(b₁, b₂, b₃, b₄ = −b₁−b₂−b₃)` and requires all six pairwise dot
products `bᵢ·bⱼ` to be ≤ 0. Its output is unique up to a 24-fold
symmetry. Used widely in crystallographic structure determination.

## What this package is for

The stated motivation of this package is twofold:

1. **DFT cell normalisation.** Electronic-structure calculations are
   more accurate and converge faster in a basis with small, nearly
   orthogonal primitive vectors. Reducing the input cell before
   passing it to a DFT code pays for itself immediately.
2. **Spacegroup analysis.** The point group of a lattice is the set
   of orthogonal transformations that map the lattice to itself. In
   a Minkowski-reduced basis, these transformations become signed
   permutation matrices (approximately, up to numerical tolerance)
   — a finite and small set to search. Without reduction first, the
   search space is unbounded. This is how
   [Spacey.jl](https://github.com/glwhart/Spacey.jl) uses this
   package.

For both use cases, what matters is that the output basis achieves
the successive minima; the absence of a canonical form (sign,
permutation, boundary ambiguities) is irrelevant.

## Further reading

- [Lattice reduction](https://en.wikipedia.org/wiki/Lattice_reduction)
  (Wikipedia) — accessible introduction with the 2D Gauss–Lagrange
  algorithm in pseudocode.
- Nguyen & Stehlé, *Low-Dimensional Lattice Basis Reduction
  Revisited*, ANTS-VI 2004 — the algorithm this package implements.
  [DOI](https://doi.org/10.1007/978-3-540-24847-7_26) ·
  [author preprint (PDF)](https://perso.ens-lyon.fr/damien.stehle/downloads/lowdim-final.pdf)
- Minkowski, *Geometrie der Zahlen*, 1910 — the original theory.
- Ryshkov & Baranovskii, *C-types of n-dimensional lattices and
  5-dimensional primitive parallelohedra* (1976) — the 3D reduction
  conditions.
- Nguyen & Stehlé, *An LLL algorithm with quadratic complexity*,
  SIAM J. Computing 2009 — the L² algorithm, relevant for
  higher-dimensional extensions.
  [DOI](https://doi.org/10.1137/070705702)
