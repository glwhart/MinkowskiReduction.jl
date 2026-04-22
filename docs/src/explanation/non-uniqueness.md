```@meta
CurrentModule = MinkowskiReduction
```

# Non-uniqueness

Minkowski reduction does **not** produce a unique reduced basis. This
is not a bug in the algorithm; it is a property of the definition.
Any algorithm that produces a Minkowski-reduced basis — this one
included — will sometimes return different-looking results for the
same physical lattice. This package does not provide canonical forms that you would see in spacegroup tables or a crystallography handbook. For those needs, see [FINDSYM](https://iso.byu.edu/findsym.php).

Three distinct sources of ambiguity coexist in 3D.

## 1. Sign ambiguity (8-fold)

Any sign flip of any basis vector produces another valid reduced
basis: if `(U, V, W)` is reduced, so is `(±U, ±V, ±W)`. That is,
eight different reduced triples define the same lattice. The
algorithm makes no deliberate choice about signs; the output signs
depend on arithmetic detail and the input.

Observed in practice:
[`DeviousMat(20)`](@ref) and [`DeviousMat(26)`](@ref) both reduce to
the simple cubic lattice, but one returns
`([0,0,−1], [1,0,0], [0,1,0])` while the other returns
`([1,0,0], [0,0,−1], [0,−1,0])` — same lattice, different signs and
axis labelling.

## 2. Permutation ambiguity at tied norms

When two or three of the successive minima are equal — which is the
normal case for *every* highly symmetric Bravais lattice (cubic,
FCC, BCC, hexagonal, tetragonal) — several orderings of the basis
vectors are all valid Minkowski-reduced bases. The outer-sort step in
`minkReduce` uses `sortperm`, whose tie-breaking is
implementation-defined, and therefore the visible column order
depends on the precise input and on the Julia version.

## 3. Boundary ambiguity at equality

If any of the twelve defining inequalities (see
[Reference → Conditions](../reference/conditions.md)) holds with
equality, the basis is *simultaneously* reducible in two different
directions: you can use the current vector, *or* you can use the
alternative that produces the equality.

For example, if `‖W‖ = ‖U + V + W‖`, then `W` and `U + V + W` are
both "shortest W for this coset" — each gives a valid reduced basis.
These equality cases mark the boundaries between regions of
configuration space where different reduced bases are canonical. On
the interior of a region (all twelve inequalities strict) the reduced
basis is unique up to sign; on a boundary, two or more choices are
simultaneously valid.

## What is unique?

Three things are intrinsic invariants of the lattice and therefore
*independent* of which Minkowski-reduced basis the algorithm happens
to return:

1. The **multiset of norms** `{‖U‖, ‖V‖, ‖W‖}` — i.e. the three
   successive minima `λ₁, λ₂, λ₃`.
2. The **orthogonality defect** `‖U‖·‖V‖·‖W‖ / |det|` — see
   [`orthogonalityDefect`](@ref).
3. The **absolute determinant** `|det|`.

These are what to compare when you want to ask "do these two bases
span the same lattice?" They are also what the package's tests
assert — see the property-based test suite in `test/runtests.jl` for
examples.

## What is not unique but could be, with normalisation

If you want a **canonical** 3D reduced basis — a deterministic
function of the lattice rather than the input basis — you need to
apply a further reduction step on top of Minkowski. Two crystallographic
standards do exactly this:

- **Niggli reduction** adds a set of sign- and permutation-normalising
  conventions that select a unique representative from each class of
  Minkowski-reduced bases. This is the standard form in
  crystallography.
- **Selling / Delaunay reduction** uses an alternative set of
  inequalities that yield a unique representation of the lattice up
  to a 24-fold symmetry group.

This package deliberately stops at Minkowski. Callers that need
uniqueness (e.g. point-group detection; database lookup of reduced
cells) should post-process with a Niggli reducer. For downstream use
cases like spacegroup analysis — which only cares about the point
group's abstract structure — the non-uniqueness here is harmless.

## Strict-inequality region: unique up to sign

If all twelve Minkowski inequalities hold *strictly* — no ties in
norm, no equalities in any of the 2D-pair or triple-linear-combination
conditions — then the reduced basis is unique **up to the 8-fold sign
ambiguity of §1**. In other words, once you are in the interior of
the fundamental domain, permutation and boundary ambiguity both
vanish; only sign remains.

For a generic (triclinic) lattice, this is the common case. For
symmetric lattices (cubic, hexagonal, etc.), you are always on a
boundary and the other ambiguities always apply.

## Practical implications

- The algorithm is **deterministic** for a given input: running
  `minkReduce(M)` twice gives the same output. The non-uniqueness is
  between **different inputs for the same lattice**.
- Two runs on the same physical structure — say, one before and one
  after a re-labelling of atoms — may give differently-signed or
  permuted reduced bases, and code that checks pointwise equality of
  output will flag a spurious mismatch.
- The algorithm-tracked transform matrix `P` correctly relates the
  input to the output that was actually produced:
  `M * P == R` always. `P` is specific to this reduction, not to the
  lattice.
