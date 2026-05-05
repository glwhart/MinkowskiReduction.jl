```@meta
CurrentModule = MinkowskiReduction
```

# The 12 Minkowski conditions

A basis `{U, V, W}` of a 3D lattice is **Minkowski reduced** if and
only if all twelve of the following inequalities hold. These are the
conditions checked by [`is_mink_reduced`](@ref), tested in the order
listed below.

For *why* these twelve conditions are sufficient (and why none of them
are redundant), see [Explanation → Algorithm](../explanation/algorithm.md).

## Conditions

| # | Inequality             | What it rules out                      | Enforced by                                       |
|---|------------------------|----------------------------------------|---------------------------------------------------|
| 1 | `‖U‖ ≤ ‖V‖`            | Mis-ordering of U and V                | Outer sort in `mink_reduce`                        |
| 2 | `‖V‖ ≤ ‖W‖`            | Mis-ordering of V and W                | Outer sort in `mink_reduce`                        |
| 3 | `‖V‖ ≤ ‖U + V‖`        | V reducible by subtracting U           | `gauss_reduce`                                     |
| 4 | `‖V‖ ≤ ‖U − V‖`        | V reducible by adding U                | `gauss_reduce`                                     |
| 5 | `‖W‖ ≤ ‖U + W‖`        | W reducible by subtracting U           | `shorten_w_in_uvw` four-corner step                |
| 6 | `‖W‖ ≤ ‖U − W‖`        | W reducible by adding U                | `shorten_w_in_uvw` four-corner step                |
| 7 | `‖W‖ ≤ ‖V + W‖`        | W reducible by subtracting V           | `shorten_w_in_uvw` four-corner step                |
| 8 | `‖W‖ ≤ ‖V − W‖`        | W reducible by adding V                | `shorten_w_in_uvw` four-corner step                |
| 9 | `‖W‖ ≤ ‖U + V + W‖`    | W reducible by subtracting U + V       | `shorten_w_in_uvw` four-corner step                |
| 10| `‖W‖ ≤ ‖U − V + W‖`    | W reducible by V − U                   | `shorten_w_in_uvw` four-corner step                |
| 11| `‖W‖ ≤ ‖U + V − W‖`    | W reducible by adding U + V            | `shorten_w_in_uvw` four-corner step                |
| 12| `‖W‖ ≤ ‖U − V − W‖`    | W reducible by U − V                   | `shorten_w_in_uvw` four-corner step                |

## Interpretation

- **Conditions 1–2** are the norm-ordering convention. They are a
  labelling choice, not a geometric property: any permutation of a
  reduced basis gives another valid basis for the same lattice.
- **Conditions 3–4** assert that in the `U-V` plane, `V` is shorter
  than any lattice vector you could get by adding or subtracting `U`.
  This is the 2D Minkowski condition that Gauss reduction enforces.
- **Conditions 5–8** say `W` cannot be shortened by `±U` or `±V`.
- **Conditions 9–12** say `W` cannot be shortened by `±(U+V)` or
  `±(U−V)` either. Together with 5–8 these cover all eight Voronoi
  neighbours of the projection of `W` onto the `U-V` plane.

Conditions 5–12 are collectively equivalent to: *`W` is the shortest
element of its coset in `L / (ℤU + ℤV)`*, which is what
[`shorten_w_in_uvw`](@ref MinkowskiReduction.shorten_w_in_uvw) enforces
in one algorithmic step.

## Tolerance

The inequalities are strict (`≤`) but [`is_mink_reduced`](@ref) tests
them up to a floating-point tolerance that scales with the largest
norm:

```julia
tol = eps(max(norm(U), norm(V), norm(W)))
norm(V) > norm(W) + tol  &&  return false  # for example
```

A bare `eps()` would be too tight for lattices whose norms are much
larger than 1 and too loose for lattices with very small norms. See
[Explanation → Precision](../explanation/precision.md#Scale-aware-tolerance-in-is_mink_reduced).

## Non-uniqueness and equality cases

If any of the twelve inequalities holds with equality, the basis is
still Minkowski reduced, but it is **not the only reduced basis for
the lattice**. For example, if `‖W‖ = ‖U + V + W‖`, then replacing `W`
with `U + V + W` yields another equally valid reduced basis. This and
other sources of non-uniqueness are discussed in
[Explanation → Non-uniqueness](../explanation/non-uniqueness.md).
