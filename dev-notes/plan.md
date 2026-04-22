# Plan: Suggested Changes to MinkowskiReduction.jl

> **Note.** This file is a reconstruction. It was written during an
> AI-assisted session as a plan document to accompany `research.md`,
> lived on the working tree for a while, and was deleted before any
> commit because its three recommendations had by then all been
> applied to the source. It is reconstructed here (from session
> context) so that the complete artefact — research analysis plus
> change proposals — is visible to future readers. All three proposals
> below were accepted and implemented in commit `d78b511` ("Docstring,
> tolerance, and error-message cleanup").

After a careful reading of the source and of the relevant literature
(see `research.md`), the core algorithm is correct and no structural
change is warranted. Three small, conservative improvements are proposed
below. All are optional.

None of these changes alters the algorithm or its output. They improve
robustness at large scale, fix a docstring that is inconsistent with the
function's return value, and correct a misleading error message.

---

## Fix 1 — Scale-aware tolerance in `isMinkReduced`

**File**: `src/MinkowskiReduction.jl`, lines 200–214.

**Problem.** Every condition is checked with an absolute tolerance of
`eps() ≈ 2.22 × 10⁻¹⁶`. The floating-point error in `norm(v)` is
proportional to `‖v‖`, so for lattices with norms ≫ 1 this tolerance is
too tight and can spuriously reject a basis that is in fact
Minkowski-reduced; for norms ≪ 1 it is too loose.

**Proposed change.** Use a tolerance that scales with the largest norm
in play. One clean way:

```julia
function isMinkReduced(U, V, W)
    tol = eps(max(norm(U), norm(V), norm(W)))
    if norm(U) > norm(V)+tol println("Condition 1 failed"); return false end
    if norm(V) > norm(W)+tol println("Condition 2 failed"); return false end
    if norm(V) > norm(U+V)+tol println("Condition 3 failed"); return false end
    # … and so on for conditions 4–12, all using `tol` in place of `eps()`
end
```

`eps(x)` returns the machine epsilon at the magnitude of `x`, which is
exactly the floating-point spacing near `norm(W)` and therefore the
right comparison tolerance for these inequalities.

**Why it matters.** The existing test suite passes because after
reduction the ordering conditions hold with a comfortable margin, not
near equality. But a user who presents a near-boundary configuration at
large scale (e.g. crystallographic cells with Å units and lattice
parameters of 10⁸ or more) could see `isMinkReduced` return `false` on
an output that `minkReduce` just produced. A scale-aware tolerance
closes that gap.

**Risk.** Low. The existing tests should continue to pass with the new
tolerance, because in every case either the inequality holds strictly
or both sides are zero-size. Worth verifying by running the full test
suite after the change.

---

## Fix 2 — Correct error message in `GaussReduce`

**File**: `src/MinkowskiReduction.jl`, line 94.

**Problem.** The message

```julia
error("GaussReduce: input vectors were linearly independent")
```

is raised when `norm(U)` becomes NaN, which happens precisely when the
input vectors are **linearly dependent** (parallel or zero). The message
says the opposite.

**Proposed change.**

```julia
error("GaussReduce: input vectors are linearly dependent")
```

**Risk.** None.

---

## Fix 3 — Repair the jldoctest in `minkReduce`'s docstring

**File**: `src/MinkowskiReduction.jl`, lines 10–13.

**Problem.** The docstring shows

```julia
julia> U = [1, 2, 3]; V = [-1, 2, 3]; W = [3, 0, 4]; minkReduce(U,V,W)
([-2.0, 0.0, 0.0], [0.0, -2.0, 1.0], [-1.0, 2.0, 3.0])
```

but `minkReduce(U,V,W)` returns a 4-tuple — the reduced basis **plus**
the iteration count. The actual output is

```
([-2.0, 0.0, 0.0], [0.0, -2.0, 1.0], [-1.0, 2.0, 3.0], 2)
```

Documenter.jl's doctest runner will flag this mismatch. (The test in
`runtests.jl` correctly checks all four values, so this is a docstring
issue only.)

**Proposed change.** Update the expected output in the jldoctest to
include the trailing `2`:

```julia
julia> U = [1, 2, 3]; V = [-1, 2, 3]; W = [3, 0, 4]; minkReduce(U,V,W)
([-2.0, 0.0, 0.0], [0.0, -2.0, 1.0], [-1.0, 2.0, 3.0], 2)
```

**Risk.** None.

---

## What is deliberately *not* being changed

- **Non-uniqueness.** Sign, permutation, and boundary ambiguity are
  intrinsic to Minkowski reduction. A user needing a canonical form
  should post-process with Niggli reduction. Adding such normalization
  here would be a new feature, not a fix, and would change existing
  outputs and tests. See `research.md` §4.
- **The 15-iteration cap in `minkReduce`.** Tight but sufficient;
  `DeviousMat(26)` (the worst representable integer stress case) hits
  exactly 15. Raising it would only hide genuine infinite-loop bugs.
- **The `floor`-based projection in `shortenW_in_UVW`.** Correct and
  consistent with the four-corner search; any replacement (e.g. `round`
  with a different corner set) would need re-derivation.
- **The `≈` stop in `GaussReduce`.** Load-bearing: without it the
  algorithm can oscillate when `‖U‖ = ‖V‖` exactly.

---

## Verification

After applying any of these fixes, confirm:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

still passes, and (for Fix 1) that `isMinkReduced` continues to return
`true` on all the high-aspect-ratio and noisy-lattice tests in
`runtests.jl`.
