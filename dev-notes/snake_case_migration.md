# snake_case migration plan

This document enumerates every identifier in `MinkowskiReduction.jl` that
would need to change to convert the public API (and internal helpers)
from camelCase / PascalCase-style names to snake_case, and lists every
file that touches each one.

**This is a breaking change for downstream users** (notably Spacey, plus
anything else that imports the package). See the "Deprecation strategy"
section at the bottom for a path that minimises pain.

---

## What changes

### Exported functions (9)

| Old name              | New name                | Defined at                     |
| --------------------- | ----------------------- | ------------------------------ |
| `minkReduce`          | `mink_reduce`           | `src/MinkowskiReduction.jl:35,99` |
| `GaussReduce`         | `gauss_reduce`          | `src/MinkowskiReduction.jl:203` |
| `RandUnimodMat2`      | `rand_unimod_mat2`      | `src/MinkowskiReduction.jl:282` |
| `RandLowerTri`        | `rand_lower_tri`        | `src/MinkowskiReduction.jl:306` |
| `DeviousMat`          | `devious_mat`           | `src/MinkowskiReduction.jl:353` |
| `isMinkReduced`       | `is_mink_reduced`       | `src/MinkowskiReduction.jl:369,396` |
| `orthogonalityDefect` | `orthogonality_defect`  | `src/MinkowskiReduction.jl:243,264` |
| `RandUnimodMat3`      | `rand_unimod_mat3`      | `src/MinkowskiReduction.jl:430` |
| `isPermutationMatrix` | `is_permutation_matrix` | `src/MinkowskiReduction.jl:488` |

### Internal (non-exported) functions (2)

| Old name           | New name             | Defined at                       |
| ------------------ | -------------------- | -------------------------------- |
| `shortenW_in_UVW`  | `shorten_w_in_uvw`   | `src/MinkowskiReduction.jl:133`  |
| `FibonacciMat`     | `fibonacci_mat`      | `src/MinkowskiReduction.jl:325`  |

### What does NOT change

These already follow Julia community convention and should be left alone:

- **Module name `MinkowskiReduction`** — Julia style is PascalCase for
  modules and types. Renaming would also force the package directory,
  `Project.toml` `name`, and every `using`/`import` statement everywhere
  to change.
- **Mathematical single-letter locals** — `U`, `V`, `W`, `M`, `P`, `δP`,
  `P_G`, `U′`, `V′`, `W′`, `a`, `b`, `c`, etc. These are math conventions
  for basis vectors / transform matrices, not camelCase. Snake_case is
  about multi-word identifiers; single capitals don't have a "case" to
  convert.
- **Local variables that are already snake_case or single-token**: `tol`,
  `denom`, `ratio`, `dists`, `case`, `norms`, `mat`, `ones_in_row`,
  `zeros_in_row`, `M_transform`, `rowvals`, `colvals`, `f1/f2/f3`, etc.

---

## File-by-file change list

### `src/MinkowskiReduction.jl`

The export line:

```julia
# before
export GaussReduce, RandUnimodMat2, RandLowerTri, minkReduce, DeviousMat, isMinkReduced, orthogonalityDefect, RandUnimodMat3, isPermutationMatrix

# after
export gauss_reduce, rand_unimod_mat2, rand_lower_tri, mink_reduce, devious_mat, is_mink_reduced, orthogonality_defect, rand_unimod_mat3, is_permutation_matrix
```

Plus rename each `function` definition, and update all internal call
sites (e.g. `minkReduce` calls itself; `shortenW_in_UVW` calls
`GaussReduce`; the matrix-form `minkReduce(M)` calls the vector form;
the matrix-form `orthogonalityDefect` calls the vector form; the
matrix-form `isMinkReduced` calls the vector form).

Also update every docstring — every `[`name`](@ref)` cross-reference,
every prose mention, and every doctest's input/output transcript needs
the new name.

### `test/runtests.jl`

About **65 call sites** across the file. All of:
`minkReduce`, `GaussReduce`, `orthogonalityDefect`, `DeviousMat`,
`isMinkReduced`, `RandUnimodMat3`, `isPermutationMatrix`,
`MinkowskiReduction.RandLowerTri`,
`MinkowskiReduction.RandUnimodMat2`,
`MinkowskiReduction.FibonacciMat`,
`MinkowskiReduction.shortenW_in_UVW`.

The `MinkowskiReduction.X` qualified calls just need the `X` part
renamed (the module name stays).

### `README.md`

Mentions `minkReduce` (5×) and `DeviousMat` (1×).

### `docs/src/index.md`

`minkReduce` (3×), `GaussReduce` (1×), `orthogonalityDefect` (1×),
`isMinkReduced` (1×).

### `docs/src/tutorial.md`

`minkReduce` (~7×), `isMinkReduced` (2×). Doctest outputs need
re-checking after the rename — the prompts (`julia>`) embed the old
names.

### `docs/src/how-to.md`

`minkReduce` (~7×), `isMinkReduced` (3×), `orthogonalityDefect` (2×),
`GaussReduce` (3×), `DeviousMat` (3×), `RandUnimodMat3` (2×),
`FibonacciMat` (1×), plus the anchor link
`#Scale-aware-tolerance-in-isMinkReduced` (line 90) which should be
updated when `precision.md`'s heading changes.

### `docs/src/explanation/algorithm.md`

`minkReduce`, `GaussReduce`, `shortenW_in_UVW`, `DeviousMat` (multiple
each), and the qualified ref
`[`shortenW_in_UVW(U, V, W)`](@ref MinkowskiReduction.shortenW_in_UVW)`
on line 54 needs both the link text and the `@ref` target updated.

### `docs/src/explanation/precision.md`

`minkReduce`, `GaussReduce`, `isMinkReduced`, `shortenW_in_UVW`,
`DeviousMat` (multiple each). Two **section headings** also embed the
old names:

- `## floor-ing near integers in shortenW_in_UVW can …` (line 35)
- `## Scale-aware tolerance in isMinkReduced` (line 92)
- `## The ≈ termination in GaussReduce` (line 126)
- `## round near half-integers in GaussReduce — non-uniqueness` (line 62)

Renaming a heading changes its anchor, so the cross-link in
`how-to.md:90` needs to be updated in lockstep.

### `docs/src/explanation/non-uniqueness.md`

`minkReduce`, `DeviousMat`, `orthogonalityDefect` (a few each).

### `docs/src/explanation/context.md`

`GaussReduce`, `orthogonalityDefect`.

### `docs/src/reference/conditions.md`

Mentions in API reference text.

### `dev-notes/plan.md`, `dev-notes/research.md`, `aihelp2-0-1.md`

Historical / planning notes that mention the API. These are optional —
keeping the old names here preserves the historical record (e.g.
`aihelp2-0-1.md` is an AI-help transcript from a prior version). Decide
per-file whether they're "living docs" that should track the rename or
"snapshots" that should not.

### `CLAUDE.md`

The "Key exported functions" section lists them by name. Update once the
rename lands so future Claude sessions see the new names.

---

## Suggested execution order

1. Land the rename in `src/MinkowskiReduction.jl` (definitions + export list + docstrings + internal call sites).
2. Update `test/runtests.jl`. Run `julia --project test/runtests.jl`; all tests should pass without any other change. (This is the cleanest sanity check before touching docs.)
3. Update `README.md`.
4. Update each file in `docs/src/`. Build docs locally (`cd docs && julia --project make.jl`) and confirm no `@ref` links are broken — Documenter will flag any cross-reference that points at an old name.
5. Update `CLAUDE.md`.
6. Bump the package version in `Project.toml` to **3.0.0** — this is a public-API break, so semver demands a major version bump.
7. Tag and announce. Notify Spacey (and any other known dependents) before tagging so they can prepare a coordinated update.

---

## Deprecation strategy (recommended)

Rather than a hard break, you can ship one transitional release that
exports both names, with the old names emitting a deprecation warning.
Julia has built-in support via `Base.@deprecate`:

```julia
# at the bottom of src/MinkowskiReduction.jl, after the new definitions
@deprecate minkReduce          mink_reduce
@deprecate GaussReduce         gauss_reduce
@deprecate RandUnimodMat2      rand_unimod_mat2
@deprecate RandLowerTri        rand_lower_tri
@deprecate DeviousMat          devious_mat
@deprecate isMinkReduced       is_mink_reduced
@deprecate orthogonalityDefect orthogonality_defect
@deprecate RandUnimodMat3      rand_unimod_mat3
@deprecate isPermutationMatrix is_permutation_matrix
```

`@deprecate` automatically:

- exports the old name,
- forwards calls to the new name,
- emits a `depwarn` the first time each old name is called.

A two-step release plan:

- **v2.1.0** — adds the new snake_case names *and* the `@deprecate`
  shims. Existing user code keeps working but warns. Spacey can migrate
  on its own schedule.
- **v3.0.0** — removes the deprecation shims. Hard break.

This is what Julia Base itself does for renames, and it's by far the
least disruptive path for downstream consumers.

> I like the two step release plan. Before we pull the trigger, could you carefully check whether there are other things in this package that break Julia style or idioms? We should do all changes now. Let's make sure not to miss anything that makes this fully julianic, if you know what I mean.

---

# Additional "make it julianic" review

I went through `src/MinkowskiReduction.jl` line-by-line and grouped the
findings by severity. Each item lists the location, the issue, and the
proposed fix. Tick the ones you want to bundle into the v2.1 release.

## Tier 1 — real issues (recommended)

### 1. `println` for failure reporting in `isMinkReduced` (`src/MinkowskiReduction.jl:374-385`)

```julia
if norm(U) > norm(V)+tol println("Condition 1 failed"); return false end
```

`println` from a library function is unidiomatic — it pollutes stdout,
can't be silenced by callers, and shows up unexpectedly in tests and
notebooks. The Julia idiom is `@debug` (off by default, opt-in via
`JULIA_DEBUG`).

Combined with snake_case + early-return-via-`||`:

```julia
function is_mink_reduced(U, V, W)
    tol = 8 * eps(max(norm(U), norm(V), norm(W)))
    norm(U) ≤ norm(V) + tol         || (@debug "Condition 1 failed";  return false)
    norm(V) ≤ norm(W) + tol         || (@debug "Condition 2 failed";  return false)
    norm(V) ≤ norm(U + V) + tol     || (@debug "Condition 3 failed";  return false)
    norm(V) ≤ norm(U - V) + tol     || (@debug "Condition 4 failed";  return false)
    norm(W) ≤ norm(U + W) + tol     || (@debug "Condition 5 failed";  return false)
    norm(W) ≤ norm(U - W) + tol     || (@debug "Condition 6 failed";  return false)
    norm(W) ≤ norm(V + W) + tol     || (@debug "Condition 7 failed";  return false)
    norm(W) ≤ norm(V - W) + tol     || (@debug "Condition 8 failed";  return false)
    norm(W) ≤ norm(U + V + W) + tol || (@debug "Condition 9 failed";  return false)
    norm(W) ≤ norm(U - V + W) + tol || (@debug "Condition 10 failed"; return false)
    norm(W) ≤ norm(U + V - W) + tol || (@debug "Condition 11 failed"; return false)
    norm(W) ≤ norm(U - V - W) + tol || (@debug "Condition 12 failed"; return false)
    return true
end
```

Users who actually want to see which condition fails do
`ENV["JULIA_DEBUG"] = MinkowskiReduction` before calling.
> Implement this change

### 2. Variable shadowing + `SparseArrays` name collision in `isPermutationMatrix` (`src/MinkowskiReduction.jl:488-506`)

Two issues in one function:

- **`rowvals` shadows `SparseArrays.rowvals`**, a real exported function.
  If a user has `using SparseArrays` loaded, intent gets murky at a
  glance. Rename to `row` and `col`.
- **The `i ->` lambda inside `count` shadows the outer `for i in 1:3`.**
  Harmless because the lambda's `i` is the inner iteration variable, but
  it makes the body confusing. Rename inner to `j`.

```julia
function is_permutation_matrix(M::AbstractMatrix{<:Real}; atol = sqrt(eps()))
    size(M) == (3, 3) || return false
    for i in 1:3
        row = @view M[i, :]
        col = @view M[:, i]
        ones_in_row  = count(j -> isapprox(abs(row[j]), 1; atol=atol), 1:3)
        zeros_in_row = count(j -> isapprox(row[j], 0;     atol=atol), 1:3)
        ones_in_col  = count(j -> isapprox(abs(col[j]), 1; atol=atol), 1:3)
        zeros_in_col = count(j -> isapprox(col[j], 0;     atol=atol), 1:3)
        (ones_in_row == 1 && zeros_in_row == 2 &&
         ones_in_col == 1 && zeros_in_col == 2) || return false
    end
    return true
end
```

(Also gets a free perf win from `@view` instead of allocating slices.)
> Explain why @view is a win? Teach me about this.

**A:** `@view M[i, :]` returns a `SubArray` — a tiny struct that holds (a) a reference to the parent `M` and (b) the indices `(i, :)`. Reading `view[k]` translates to `M[i, k]` *at access time*. No memory is copied.

`M[i, :]` *without* `@view` does the opposite: it allocates a fresh `Vector{T}` and copies the row's elements into it. For a 3-row matrix that's a small heap allocation plus a copy, every time.

In `is_permutation_matrix` the for-loop runs 3 times, and each iteration without `@view` allocates two 3-element vectors (`row` and `col`). So **6 small heap allocations per call** vs **0** with `@view`. Microscopic in absolute terms (≲ 200 bytes, microsecond-ish), but:

1. It's the *idiomatic* way to slice when you only need to read.
2. It scales: if you ever generalise this to higher-dimensional lattices, those allocations grow with the dimension.
3. It removes a class of bug that doesn't apply here but is worth being aware of (see caveat below).

**Two caveats to know about views:**

- **A view aliases the parent.** If the parent is later mutated, the view reflects the new values. In `is_permutation_matrix` we only read, so this is fine. But code like `r = @view M[1, :]; M[1,1] = 0; use(r)` would surprise.
- **Some old code paths dispatch differently on `Vector` vs `SubArray`.** Mostly irrelevant in modern Julia but worth knowing if you see a perf cliff.

**`@view` vs `@views`:** `@view` modifies one slicing expression. `@views` is a block-level macro that rewrites *every* `[…]` slice inside the block to use `view`. See item 12 for that.

### 3. `Int64` vs `Int` inconsistency

`Matrix{Int}(I, 3, 3)` is used in `minkReduce` and `GaussReduce`, but
`Matrix{Int64}(I, 3, 3)` in `RandUnimodMat3` (line 432). And
`round(Int64, …)` and `convert(Int64, …)` in `DeviousMat` /
`FibonacciMat`.

Julia's `Int` is `Int64` on 64-bit and `Int32` on 32-bit; using the
alias is the portable, idiomatic choice. The only places `Int64` is
*load-bearing* are the `DeviousMat` overflow comments — and even there,
on a 32-bit platform we'd want a sharper failure mode, not silent wrap.

Proposal: replace every `Int64` in `src/` with `Int`, except the
docstring of `DeviousMat` which legitimately documents the 64-bit
overflow boundary.
> Implement suggestion

### 4. `convert(Int64, x)` → `Int(x)` (`src/MinkowskiReduction.jl:356`)

```julia
a, b = convert(Int64, (u+v+1)/2), -u
```

Modern Julia: `Int((u+v+1)/2)` (or, since `(u+v+1)` is integer and we're
dividing by 2, `(u+v+1) ÷ 2` skips the float round-trip entirely).
> I like the latter

### 5. `using Random` is unused in `src/`

`rand` is exported from `Base`. The src never calls `Random.shuffle!`,
`MersenneTwister`, or `seed!`. Drop the import; `using LinearAlgebra` is
all that's needed.
> Implement suggestion

(The test file uses `shuffle!`, so leave its `using Random` alone.)

### 6. `orthogonalityDefect(M)` only validates one dimension (`src/MinkowskiReduction.jl:265`)

```julia
size(M, 2) == 3 || error("Matrix must have exactly 3 columns")
```

A 2×3 matrix would pass this check and then crash inside `cross` with a
less helpful error. Validate the full shape:

```julia
size(M) == (3, 3) || throw(DimensionMismatch("expected a 3×3 matrix, got $(size(M))"))
```

Same gap exists in `isMinkReduced(M)` (`:396`) — it has *no* size check
at all and silently uses the first three columns.
> Fix both of these.

### 7. Use `DimensionMismatch` / `ArgumentError` / `DomainError` instead of generic `error(...)`

Five generic `error(...)` calls throw `ErrorException`. Julia idiom is
to use the most specific built-in exception type — it lets callers
`try`/`catch` precisely:

| Location | Current | Proposed |
| --- | --- | --- |
| `:63` | `error("minkReduce: Too many iterations")` | `error("…")` is fine here (genuine bug indicator) — or define a package-local `ConvergenceError`. |
| `:150` | `error("shortenW_in_UVW: U and V are linearly dependent")` | `throw(ArgumentError("U and V are linearly dependent"))` |
| `:220` | `error("GaussReduce: input vectors are linearly dependent")` | `throw(ArgumentError(...))` |
| `:226` | `error("GaussReduce: Too many iterations")` | as `:63` |
| `:265` | `error("Matrix must have exactly 3 columns")` | `throw(DimensionMismatch(...))` |
| `:329` | `error("Overflow in FibonacciMat function")` | `throw(OverflowError("FibonacciMat: k=$k overflows Int64"))` |
| `:354` | `error("for DeviousMat, n > 2")` | `throw(DomainError(n, "DeviousMat requires n ≥ 3"))` |
| `:431` | `error("RandUnimodMat3: k must be positive")` | `throw(DomainError(k, "k must be ≥ 1"))` |

⚠️ **Test impact:** `runtests.jl` has several `@test_throws ErrorException
…` lines (e.g. `:76`, `:77`, `:79`, `:80`, `:217`). Those will need
their expected exception types updated, or relaxed to
`@test_throws Exception`.
> Implement. Thanks.

## Tier 2 — style/consistency cleanup (recommended but minor)

### 8. Mixed `∈` vs `in` in `for` loops

- `RandUnimodMat2`: `for i ∈ 1:n` (line 284)
- `RandUnimodMat3`: `for _ in 1:k` (line 434)
- `isPermutationMatrix`: `for i in 1:3` (line 493)

Pick one and apply consistently. Both are valid; `in` is more common in
Base. Recommend standardising on `in`.
> Standardize on ∈

### 9. Single-line `if cond expr end` patterns

```julia
if norm(U) > norm(V) || norm(U) ≈ norm(V) break end          # line 225
```

Idiomatic Julia uses short-circuit `&&` for one-liners:

```julia
(norm(U) > norm(V) || norm(U) ≈ norm(V)) && break
```

The 12 `if … println … return false end` lines in `isMinkReduced` get
fixed by item 1 above.
> I like your suggestion

### 10. `m == 0 && (m = 1)` (`:441`, `:451`)

Works, but the parenthesised assignment as a side-effect of `&&` is
mildly cryptic. Cleaner:

```julia
m = iszero(m) ? 1 : m
```

Trivial. Take it or leave it.
> Take it.

### 11. `prod(norm.([a,b,c]))` allocates a 3-element vector (`:244`)

Idiomatic non-allocating form:

```julia
prod(norm, (a, b, c))
```

(`prod` accepts a function as first arg, just like `sum`.) Same numeric
result, no temporary array. Picosecond-level perf win, but tidier.
> Implement

### 12. `M[:,1]`, `M[:,2]`, `M[:,3]` allocate in matrix-form wrappers (`:100`, `:266`, `:397`)

Use `@views` to avoid the three-column copy in:

```julia
function mink_reduce(M)
    @views U, V, W, P, _ = mink_reduce(M[:,1], M[:,2], M[:,3])
    return hcat(U, V, W), P
end
```

(`float(view)` still does the right thing inside the inner method.)
Trivial perf win, also tidier.
> Tell me more about why `@views` is better

**A:** Same mechanism as `@view` (see item 2), different scope. `@views` is a *block-level* macro: it rewrites **every** slicing expression inside the block, function, or `let` it wraps to use `view` instead of copy.

```julia
@views function mink_reduce(M)
    U, V, W, P, _ = mink_reduce(M[:,1], M[:,2], M[:,3])  # all three become view(M, :, k)
    return hcat(U, V, W), P
end
```

Without `@views`, each `M[:,k]` allocates a fresh 3-element `Vector` and copies the column. With `@views`, all three become zero-copy `SubArray`s. For a 3×3 matrix that's three avoided allocations per call — small but free.

**Why `@views` over sprinkling `@view` everywhere:**

- **One macro at the function top, no clutter at each call site.** Easier to read, harder to miss one.
- **Consistency.** If a future contributor adds another slice expression inside the function, they get the view treatment automatically.

**When *not* to use `@views`:** if some slice expression in the function is *meant* to produce an independent copy (e.g. you mutate the slice and want the parent untouched), `@views` will silently change the semantics. Easy to spot — there are no such cases in `mink_reduce(M)`, `orthogonality_defect(M)`, or `is_mink_reduced(M)`; they all just hand the slices off to a method that reads them.

**Honest perf assessment for this package:** for 3×3 matrices the speedup is unmeasurable. The real reasons to take it are (a) it's the idiomatic style, and (b) the package's own docs already pitch it as a possible high-dimensional generalisation, and at higher dimensions the allocations matter.

### 13. Docstring formatting — `(U,V,W)` → `(U, V, W)`

Several docstring signatures have no spaces after commas:
`isMinkReduced(U,V,W)` (line 362), `orthogonalityDefect([1,1,0],…)`
(line 239), etc. Trivial; harmonise with the rest.
> Leave these

## Tier 3 — judgment calls (your decision)

### 14. Should `RandLowerTri` be exported?

Currently the export list contains `RandLowerTri` (a 2-line helper that
returns `[1 0; rand(-n:n) 1]`) but **not** `FibonacciMat` (which is more
useful for stress testing). Both are similar building blocks for
`RandUnimodMat2`.

Three coherent choices:

- **A.** Un-export `RandLowerTri` (likely nobody calls it externally).
> do this

- **B.** Also export `FibonacciMat` (consistency the other way).
- **C.** Leave both as they are.

I'd vote **A** — `rand_lower_tri` doesn't pull its weight as part of the
public surface. (If you take this, the deprecation shim from item §
"Deprecation strategy" should still cover `RandLowerTri` to avoid
breaking anyone who *was* using it.)

### 15. Type annotation consistency

Only `RandUnimodMat3(k::Integer = 10)` constrains its argument. The
others (`DeviousMat(n)`, `FibonacciMat(k)`, `RandUnimodMat2(n)`,
`RandLowerTri(n)`) take untyped `n`/`k`. For these all-integer
generators, adding `::Integer` is mildly nice for error messages but
not load-bearing.

I'd lean toward leaving the no-annotation versions alone — Julia idiom
is "don't over-constrain; let duck typing work."
> Agreed

### 16. Greek/prime unicode in identifiers (`δP`, `U′`, `V′`, `W′`)

Pure stylistic call — these are valid Julia and read naturally as math.
The Julia community is split: numerical packages (e.g. `DifferentialEquations.jl`) embrace unicode; Base avoids it.

Keep them. They make the code closer to the math notation in the
Nguyen–Stehlé paper.
> Keep them

### 17. `error("minkReduce: Too many iterations")` includes the function name

This is a deliberate convention from the original code, useful for stack-traceless errors. Either prefix-everything or prefix-nothing — currently the package does both. I'd standardise on prefixing (the `error(funcname: msg)` pattern), since it survives `catch err` chains.
> Prefix all

## Things I checked and don't think need changing

- **Module name `MinkowskiReduction`** — Julia convention is PascalCase
  for modules, so this is already right.
- **Single-letter math identifiers** (`U`, `V`, `W`, `M`, `P`) — not
  camelCase; just math notation. Leave alone.
- **Unicode operators `⋅`, `×`, `≈`, `√`** — idiomatic in Julia
  numerical code. Leave alone.
- **`while true … break`** loops — unusual in some languages but
  perfectly idiomatic Julia. Leave alone.
- **Docstring jldoctests** — well-structured, run as part of `Pkg.test`.
  Just need the renamed-identifier sweep.
- **Project layout** (`src/Module.jl`, `test/runtests.jl`,
  `Project.toml`) — standard.
- **`@inbounds` in `RandUnimodMat3`** — fine as-is; the indices are
  proven in 1:3 by construction.

## Summary checklist for v2.1 — confirmed

Reflecting your inline decisions above. Everything in this list is approved unless noted.

- [ ] snake_case rename + `@deprecate` shims (the original migration plan)
- [ ] `println` → `@debug` in `is_mink_reduced` (item 1) ✅
- [ ] `is_permutation_matrix`: rename `rowvals`/`colvals` → `row`/`col`, fix `i`-shadowing → `j`, use `@view` (item 2) ✅
- [ ] `Int64` → `Int` throughout `src/` except `DeviousMat`'s docstring overflow note (item 3) ✅
- [ ] `convert(Int64, (u+v+1)/2)` → `(u+v+1) ÷ 2` (item 4) ✅
- [ ] Drop `using Random` from `src/` (item 5) ✅
- [ ] Validate full `(3,3)` shape in `orthogonality_defect(M)` *and* `is_mink_reduced(M)` (item 6) ✅
- [ ] Specific exception types (`ArgumentError`, `DomainError`, `DimensionMismatch`, `OverflowError`) — **and update `@test_throws` lines in tests** (item 7) ✅
- [ ] Standardise on `∈` for `for`-loops in `src/` (item 8) ✅ *(your call — going with `∈` not `in`)*
- [ ] Replace single-line `if … end` with short-circuit `&&` / `||` (item 9) ✅
- [ ] `m == 0 && (m = 1)` → `m = iszero(m) ? 1 : m` in `rand_unimod_mat3` (item 10) ✅
- [ ] `prod(norm.([a,b,c]))` → `prod(norm, (a,b,c))` (item 11) ✅
- [ ] `@views` on the three matrix-form wrappers (item 12) ✅
- [ ] Docstring comma spacing — **leave alone** (item 13) ❌
- [ ] Un-export `RandLowerTri` (option **A** in item 14); deprecation shim still covers it ✅
- [ ] Type annotations — **leave alone** (item 15) ❌ *(no change)*
- [ ] Unicode `δP`, `U′`, `V′`, `W′` — **keep** (item 16) ❌ *(no change)*
- [ ] Standardise on prefixed errors (`"funcname: msg"`) (item 17) ✅

## Open questions before I start

1. **PR shape.** Single PR or two (rename first, idiom cleanup second)?
   - **Single PR pro:** one round of test-fixing, one CI run, one tag.
   - **Two PRs pro:** the rename diff stays mechanical and easy to review; the idiom changes get isolated diffs you can scrutinise. Easier to bisect later if something breaks.
   - **My recommendation:** two PRs. First PR is pure mechanical rename + deprecation shims (touches lots of files but every change is the same kind). Second PR is the idiom cleanup (smaller, but each change is a real semantic edit worth eyeballing).
I'll follow your recommendation

2. **`@deprecate` for the un-exported `RandLowerTri` (item 14).** `Base.@deprecate` exports the old name by default. We can pass `false` as a third argument to keep it un-exported: `@deprecate RandLowerTri rand_lower_tri false`. Confirm you want that, since the current state has `RandLowerTri` *exported*.
Confirmed

3. **`dev-notes/plan.md`, `dev-notes/research.md`, `aihelp2-0-1.md`** — update these (living docs) or freeze them with old names (historical record)? My suggestion: freeze `aihelp2-0-1.md` (it's a transcript), update `plan.md` / `research.md` if they're being actively edited.
freeze

4. **Version number.** v2.1.0 (additive — old names still work via deprecation) is correct semver for the first PR/release. v3.0.0 for the eventual shim removal.
Agreed

Give me the go-ahead on PR shape (and #2/#3), and I'll start implementing.

