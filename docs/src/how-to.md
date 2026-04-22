```@meta
CurrentModule = MinkowskiReduction
DocTestSetup = quote
    using MinkowskiReduction
    using LinearAlgebra
end
```

# How-to guides

Short, self-contained recipes for common tasks. Each assumes you have
already done the [tutorial](tutorial.md) and know the basic shape of
`minkReduce`'s return value. Pick the recipe you need; there is no
order.

---

## How do I get the integer transform matrix, not just the reduced basis?

The matrix form of `minkReduce` returns both:

```jldoctest howto
julia> M = [1.0 0 1; 0 1 1; 0 0 1];

julia> R, P = minkReduce(M);

julia> M * P == R
true
```

`P` is the exact integer (`Matrix{Int}`), with `|det(P)| = 1`, that converts `M` into its reduced form `R`. See
[Explanation → Algorithm](explanation/algorithm.md) for how it is
constructed.

---

## How do I transform atomic positions when I reduce the cell?

If your atoms are expressed in **fractional coordinates** relative to
the old basis `M`, their fractional coordinates in the reduced basis
`R` are `x_new = P⁻¹ · x_old`. (Cartesian coordinates are unchanged by
a change of basis.)

A minimal example with two atoms at `(0.1, 0.2, 0.3)` and `(0.5, 0.5,
0.5)` in fractional coordinates:

```jldoctest howto
julia> M = [1.0 0 1; 0 1 1; 0 0 1];

julia> X_old = [0.1 0.5; 0.2 0.5; 0.3 0.5]        # columns = atomic positions
3×2 Matrix{Float64}:
 0.1  0.5
 0.2  0.5
 0.3  0.5

julia> R, P = minkReduce(M);

julia> X_new = mod.(inv(P) * X_old, 1)            # transform, then wrap into [0, 1)
3×2 Matrix{Float64}:
 0.4  0.0
 0.5  0.0
 0.3  0.5
```

The `mod.(..., 1)` wrap is the standard crystallographic convention:
fractional coordinates of `1.0` label the same lattice site as `0.0`,
so both atom 2's transformed components of `1.0` collapse to `0.0`.
`X_old` and `X_new` describe the same physical atomic positions,
expressed in different bases.

`inv(P)` is exact for unimodular matrices — its entries are always
integers, recoverable as typed integers via `round.(Int, inv(P))` if
you need them that way.

---

## How do I check whether a basis is already Minkowski reduced?

```jldoctest howto
julia> isMinkReduced([1.0 0 0; 0 1 0; 0 0 1])
true
```

There is also a three-vector form `isMinkReduced(U, V, W)`. On a
non-reduced input the function returns `false` and also prints which of
the 12 conditions failed (see
[The 12 Minkowski conditions](reference/conditions.md)).

The test uses a scale-aware floating-point tolerance — see
[Explanation → Precision](explanation/precision.md#Scale-aware-tolerance-in-isMinkReduced).

---

## How do I compute the orthogonality defect?

```jldoctest howto
julia> orthogonalityDefect([1.0 0 0; 0 1 0; 0 0 1])   # orthogonal
1.0

julia> orthogonalityDefect([1.0 1 0; 1 0 1; 0 1 1])   # FCC conventional basis
1.4142135623730954
```

The defect is the ratio `‖a‖·‖b‖·‖c‖ / |det|`. It equals 1 exactly for
an orthogonal basis and grows without bound as the basis becomes more
skewed. A Minkowski-reduced basis minimizes this quantity over all
bases of the lattice.

---

## How do I reduce a 2D sublattice?

Use `GaussReduce` directly. It operates on any pair of vectors — the
vectors can themselves live in 3D, which is useful for reducing an
in-plane sublattice of a crystal (e.g. for a slab calculation):

```jldoctest howto
julia> U = [5.0, 8.0, 0.0]; V = [8.0, 13.0, 0.0];

julia> a, b, P = GaussReduce(U, V);

julia> (a, b)
([0.0, -1.0, 0.0], [-1.0, 0.0, 0.0])

julia> hcat(U, V) * P == hcat(a, b)
true
```

The return is `(shorter, longer, 2×2 integer transform)`. See
[Explanation → Algorithm](explanation/algorithm.md) for why Gauss
reduction is the engine behind `minkReduce`.

---

## How do I verify `P` is unimodular exactly?

For a small `P` with entries of magnitude ~1, `abs(det(P)) ≈ 1` via
Float64 is exact. This will be the case for most materials science problems. For pathological `P` (e.g. when reducing a heavily skewed
cell — see [`DeviousMat`](@ref)), Float64's `det` loses precision and
can return a number far from `1.0` even though the true value is
exactly `1`. Use `BigInt` for an exact check:

```jldoctest howto
julia> M = DeviousMat(26);

julia> R, P = minkReduce(M);

julia> abs(det(BigInt.(P))) == 1
true
```
This should never be needed for typical cases. Pathological cases like this are helpful to make the algorithm robust.
See [Explanation → Precision](explanation/precision.md#det-P-loses-precision-in-Float64)
for why this matters.

---

## How do I generate test inputs for the reducer?

Three utilities ship with the package:

- [`RandUnimodMat3(k)`](@ref) — a random `3×3` integer matrix with
  `|det| = 1`. Larger `k` produces larger entries.
- [`DeviousMat(n)`](@ref) — a heavily-disguised simple-cubic basis,
  the worst-case 3D input for the reducer. Bounded by `3 ≤ n ≤ 26`;
  larger `n` silently overflows `Int64`.
- `FibonacciMat(k)` (unexported) — an ill-conditioned 2D matrix with
  consecutive Fibonacci numbers as entries. Overflows around
  `k ≈ 92`. Used in some unit tests.

```jldoctest howto
julia> M = RandUnimodMat3(5);

julia> size(M)
(3, 3)

julia> abs(det(BigInt.(M))) == 1
true

julia> size(DeviousMat(26))
(3, 3)
```

---

## How do I handle linearly-dependent input?

`minkReduce` and `GaussReduce` detect linear dependence by watching for
`NaN` norms during the iteration and raise `ErrorException`:

```jldoctest howto
julia> U = [2.0, 0, 0]; V = [0.0, 2, 0]; W = U + V;  # W lies in U-V plane

julia> try
           minkReduce(U, V, W)
       catch e
           e.msg
       end
"GaussReduce: input vectors are linearly dependent"
```

In practice this fires when the input is mathematically degenerate **or
numerically indistinguishable from degenerate** — e.g. when `W` lies
within ~10⁻¹⁷⁰ of the U-V plane (for unit-scale vectors). Slightly
larger perturbations like `1e-150` are handled correctly.

If this error triggers, it most likely means the input is bad (vectors are linearly dependent and don't define a basis), not that finite precision or a bug caused the error. The user should check the input.
