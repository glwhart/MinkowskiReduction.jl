```@meta
CurrentModule = MinkowskiReduction
```

# Floating-point behaviour

Low-dimensional Minkowski reduction does **not** suffer the dramatic
cancellation pathologies of high-dimensional LLL — in 2D and 3D the
Gram-Schmidt coefficients are bounded and the algebra is shallow.
Nevertheless several details of finite-precision behaviour are worth
understanding if you are debugging unexpected output, feeding the
reducer extreme input, or extending it.

## Inputs are promoted to `Float64` on entry

`mink_reduce` and `gauss_reduce` both begin by promoting their input
vectors to `Float64`:

```julia
U, V, W = float(U), float(V), float(W)
```

This exists to avoid Int64 overflow in intermediate dot products. For
inputs like [`devious_mat(26)`](@ref) (used for stringent testing), whose entries are ~10¹⁴, a naïve
integer dot product `U⋅V` would overflow Int64 (max ≈ 9.2 × 10¹⁸)
well before the algorithm's first rounding decision. The `Float64`
representation is exact for integers up to `2⁵³ ≈ 9 × 10¹⁵`, so
`devious_mat(26)` is representable without loss; larger inputs lose
mantissa bits but do not throw exceptions.

The transform matrix `P` is tracked in parallel as `Matrix{Int}` and
is *not* affected by the float promotion — every operation applied to
the vectors is mirrored as an exact integer column operation on `P`.

## `floor`ing near integers in `shorten_w_in_uvw` can increase iteration count but algorithm will still succeed

The projection coefficients

```julia
a = floor(Int, ((U⋅W)(V⋅V) − (V⋅W)(U⋅V)) / denom)
b = floor(Int, ((V⋅W)(U⋅U) − (U⋅W)(U⋅V)) / denom)
```

can be wrong by 1 when the argument is a near-integer floating-point
value. If the exact value of the projection is `3.0` but the computed
value is `2.9999999…` due to roundoff, `floor` returns `2` rather
than `3`, and `W` is translated into an adjacent parallelogram. The
four-corner search then finds a corner that is shorter than the
input but not the global optimum for this iteration.

The outer loop corrects this on the next pass, because the sum
`‖U‖² + ‖V‖² + ‖W‖²` still decreases strictly. Correctness is
preserved; only the iteration count is affected. This is why the
iteration cap exists for an algorithm that provably terminates.
Empirical cost of one such off-by-one is ≤ 1 extra iteration per
occurrence. `devious_mat(26)` — the worst integer input —
converges in exactly 15 iterations; the cap is set at 29, leaving
room for a modest accumulation of Float64 corrections on top.
See [Algorithm → Layer 3](algorithm.md#layer-3-the-outer-loop)
for the derivation.

## `round` near half-integers in `gauss_reduce` — non-uniqueness

Inside `gauss_reduce` the update is

```julia
m = round(Int, (U⋅V) / (U⋅U))
V, U = U, V − m·U
```

When `(U⋅V)/(U⋅U) ≈ ±0.5`, Julia's banker's rounding (round-half-to-
even) can go either way, producing **two different but equally valid
reduced pairs**. This is the floating-point manifestation of the
equality case in conditions 3 and 4 — a genuine non-uniqueness of
Minkowski reduction, not a bug. See [Non-uniqueness](non-uniqueness.md).

## `denom → 0` when `U ∥ V`

```julia
denom = (U⋅U)(V⋅V) − (U⋅V)²
```

is the squared area of the (U, V) parallelogram; it vanishes when
`U` and `V` are parallel. Floating-point division then produces `Inf`
or `NaN`, which propagates into the `floor` arguments. The algorithm
detects this via an `isfinite` check on the ratios and raises
`ErrorException("… linearly dependent")`. The failure is deterministic
at the scale at which numerical linear dependence begins, which in
practice is around `‖cross‖ / (‖U‖·‖V‖) ≲ 10⁻¹⁷⁰` for unit-scale
inputs. Perturbations of `1e-150` and above are handled without error.

## Scale-aware tolerance in `is_mink_reduced`

Each of the 12 conditions is tested as

```julia
tol = eps(max(norm(U), norm(V), norm(W)))
if norm(U) > norm(V) + tol ... end
```

The tolerance `eps(x)` is the spacing of `Float64` numbers near the
magnitude of `x` — i.e. the natural yardstick for comparing two
quantities of that size. A bare `eps()` (≈ 2.2 × 10⁻¹⁶) would be the
right tolerance only if all norms are ~1; it would be too tight for
a lattice with norms of ~10⁸ (where floating-point error in
`norm(v)` is ~10⁻⁸ and a strict `>` test spuriously fails) and too
loose for a lattice with norms of ~10⁻⁸.

The test suite exercises aspect ratios from 10⁻⁸ to 10⁸ without
spurious `is_mink_reduced` failures thanks to this scaling.

## `det(P)` loses precision in Float64

Even though `P` is exact integer and `|det(P)| = 1` in every case,
Julia's `det(::Matrix{Int})` internally promotes to `Float64` via LU
factorisation and can return a number far from 1 for large `P`. For
example, the `P` from `mink_reduce(devious_mat(26))` has entries of
magnitude ~10¹⁴, and `det(P)` evaluates to ~10¹² — not because `P` is
near-singular, but because the 3×3 determinant involves subtractions
of terms ~10²⁸ computed in Float64.

`det(BigInt.(P))` is exact and always returns `±1`. The package's
test suite uses this pattern and the API reference's how-to guide
points callers at it.

## The `≈` termination in `gauss_reduce`

Termination of the Gauss inner loop is governed by:

```julia
if norm(U) > norm(V) || norm(U) ≈ norm(V) break end
```

The `≈` clause looks like a redundant tolerance check but is essential
for termination on symmetric inputs. See
[Algorithm](algorithm.md#Layer-1-2D-Gauss-reduction) for the
complete oscillation argument. `Base.isapprox`'s default `rtol =
sqrt(eps)` ≈ 1.5 × 10⁻⁸ is well-suited to this purpose: in 2D the
norm ratios near termination are either ~1 (the case we need to
accept) or bounded away from 1 by far more than `√eps`.

## What is *not* a problem

- **Cancellation in the algebraic projection.** The coefficient
  formulas are polynomial in dot products; no subtraction of nearly
  equal large quantities occurs in 3D.
- **Integer inputs.** For integer input (within `Int64` range), all
  dot products are mathematically exact. Float64 promotion happens
  for storage but the values are exact integers.
- **Scale invariance of iteration count.** Uniformly scaling the
  basis does not change how many iterations the algorithm takes
  (verified by the property-based tests). All algorithmic decisions
  depend only on ratios, which are scale-free.
