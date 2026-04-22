# AI-assisted session log — MinkowskiReduction.jl, work toward v2.0.1

A chronological summary of the prompts given to the assistant and the
work produced in response, covering the v1.5.0 → v2.0.1 arc.

---

## 1. Codebase orientation

**You**: `/init` — create a CLAUDE.md.
**Me**: Explored the repo and wrote `CLAUDE.md` with build/test commands and
a concise architecture overview (single-file implementation, three algorithmic
layers, key exported functions, test-suite structure).

**You**: Read the code carefully. Research Minkowski reduction, the
Nguyen–Stehlé paper, non-uniqueness, and finite-precision issues. Write
`research.md` and, if warranted, `plan.md`.
**Me**: Read the source and tests end-to-end; researched the Divio/Diátaxis
framework of documentation and the Nguyen–Stehlé 2004 paper via web. Produced:
- `research.md` — algorithm walkthrough (GaussReduce → shortenW_in_UVW →
  minkReduce), the 12 Minkowski conditions, three kinds of non-uniqueness
  (sign, permutation, boundary), and finite-precision issues (scale-blind
  tolerance, `floor` near integers, `round` at half-integers, near-parallel
  degeneracy).
- `plan.md` — three small, conservative fixes recommended: scale-aware
  tolerance in `isMinkReduced`, corrected "linearly independent"→"dependent"
  error message in `GaussReduce`, and a stale jldoctest for `minkReduce`.

**You**: What do you mean by "load-bearing"?
**Me**: Explained: the `≈` termination branch in `GaussReduce` looks like
a redundant tolerance check but is essential — without it, `‖U‖ = ‖V‖`
inputs (common in FCC/BCC/hex) cause an infinite oscillation. Asked to be
more verbose in the comments; expanded both mentions of "load-bearing" in
`research.md` to spell out the oscillation mechanism and cite concrete
lattice families where it bites.

---

## 2. Apply the three recommended fixes

**You**: Make the three suggested changes. Update docstrings and comments
accordingly.
**Me**: Applied all three: scale-aware `tol = eps(max(norm(U),norm(V),norm(W)))`
in `isMinkReduced`; typo fix in `GaussReduce`'s NaN branch; corrected the
`minkReduce` jldoctest to show the 4-tuple return. All 64,377 tests still
passed.

**You**: Are there missing or too-terse docstrings?
**Me**: Surveyed every exported and internal function. Found: no docstring
for `isMinkReduced(M)`; terse/incorrect ones for `minkReduce(M)`,
`shortenW_in_UVW`, `GaussReduce`, `RandUnimodMat2`, `RandLowerTri` (which
incorrectly described the matrix form), `FibonacciMat`, `DeviousMat`; and
broken `@ref` links to the non-existent `RandUnimodMat`. Fixed all of them.

---

## 3. Property-based test suite

**You**: Are there unit tests that could be added to make the code more
robust?
**Me**: Proposed five property-based test classes: `minkReduce` invariants
on random unimodular transforms; idempotence; `GaussReduce`-specific tests
(it had none); sign/permutation invariance of the sorted norm multiset;
scale invariance of the iteration count.

**You**: Add the first 5.
**Me**: Added them; found two bugs in the process — my docstring update had
the `GaussReduce` return order backwards (first entry is *shorter*, not
longer), and the `|det|` check needed a looser `rtol` because cancellation
in `det` loses precision on large-entry integer inputs. Test count grew
64,377 → 65,443.

**You**: How would you implement test #6 (Bravais lattices)?
**Me**: Described a parametric table approach (seven Bravais primitive
cells with precomputed sorted-norm and orthogonality-defect invariants,
each stress-tested through random unimodular transforms), with helpers
for rhombohedral cells.

**You**: Yes, add it.
**Me**: Added test #6. Count 65,443 → 65,877, all passing.

---

## 4. Spacey integration: considerations + transform matrix P

**You**: This code will be used to make Spacey more robust. Any
considerations before moving on?
**Me**: Seven points, ranked by impact. The top one: `minkReduce` should
return the integer unimodular matrix P that records the change of basis
— downstream code needs it to transform atomic positions, symmetry
operations, etc. Also: decide canonicalisation policy; keep tolerances
consistent; handle handedness; think about near-Bravais-boundary stability.

**You**: Add P to the output. Spacey doesn't need a canonical basis.
**Me**: Threaded the transform through all three layers: `GaussReduce →
(a, b, P)`; `shortenW_in_UVW → (U, V, W, δP)`; `minkReduce(U,V,W) → (U, V,
W, P, i)`; `minkReduce(M) → (R, P)`. Each column operation on the vectors
is mirrored as a right-multiplication on P, so P stays exact-integer.
Along the way:
- discovered a silent `Int64` overflow pathology for `DeviousMat(26)`-sized
  inputs. The old code had accidentally worked because mixed-type arithmetic
  promoted to `Float64`; my cleanup with `round(Int, …)` broke the escape
  hatch. Fixed by `float(U), float(V), float(W)` at function entry.
- `round(Int, NaN)` throws `InexactError` instead of propagating and
  tripping the "linearly dependent" error. Added `isfinite` guards.
- `det(P)` for `DeviousMat(26)` returns `~10¹²`, not 1, because Julia's
  `det` on `Matrix{Int}` promotes to Float64 and loses precision. The P is
  actually exact-integer unimodular — `abs(det(BigInt.(P))) == 1` confirms.
  Tests updated accordingly.

Test count grew 65,877 → 67,406. Added dedicated "test #7" for P contract.

---

## 5. v2.0.0 release

**You**: Push changes and tag a new release.
**Me**: Discussed version bump (breaking change → 2.0.0 per SemVer),
plan.md deletion, CLAUDE.md commit, commit structure (three logical
commits), and drafted release notes. You approved; I executed three
commits, pushed, tagged v2.0.0, and pushed the tag. `gh` CLI was not
installed so the GitHub Release page required manual creation.

**You**: Watch the first CI run and report.
**Me**: Polled the GitHub Actions API unauthenticated; all workflows
green on `079f839`.

---

## 6. CI for doctests and documentation

**You**: Help set up CI to exercise jldoctests and the documentation build.
**Me**: Diagnosed: `.github/workflows/Documentor.yml` was empty (0 bytes)
— docs CI had silently never been running. Found four stale jldoctests
(four doctest failures when I enabled the build locally). Proposed:
populate the workflow, fix the four doctests, delete the legacy
`CI.yml.old`, delete outdated `docs/Manifest.toml` + gitignore it, and
add `doctest(MinkowskiReduction)` to `runtests.jl` so `Pkg.test()`
validates examples too.

**You**: Proceed with the plan, doctest in runtests.
**Me**: Executed all six steps. All tests still passing, doctests now
exercised by both the docs workflow and `Pkg.test()`. Committed, pushed;
CI green.

---

## 7. Documentation expansion — Diátaxis

**You**: Read the Divio/Diátaxis documentation framework and plan a
significant expansion of the docs for MinkowskiReduction.jl.
**Me**: Fetched the Divio site and Diátaxis summary. Reported my
understanding: four quadrants (tutorial / how-to / reference /
explanation), the two axes that define them, the don't-mix-quadrants
rule, the "start with what's missing" heuristic. Wrote `docs-plan.md`
with a full target structure — 9 pages across the four quadrants — and
a per-page content outline.

**You**: Execute the plan.
**Me**: Wrote all 9 pages (tutorial, how-to, reference × 2, explanation
× 4, landing page rewrite), updated `docs/make.jl`'s nav tree, migrated
`research.md` content into the four Explanation pages, deleted the
now-duplicated `research.md`. About 1,100 lines of new markdown. Fixed
four doctest failures exposed by the new pages (iteration count, atomic-
position transform output, error message wording, random-matrix filter
syntax). Committed as `8cdcb54`, pushed; CI green.

---

## 8. Visual aids + external references

**You**: Incorporate links to the Wikipedia lattice-reduction article
and its images where helpful. Maybe add the Nguyen–Stehlé PDF to the
repo.
**Me**: Fetched both Wikipedia and Wikimedia Commons. The article has
one illustrative SVG (`Lattice-reduction.svg`), **public domain** —
better than CC-BY-SA. On the PDF: flagged copyright concerns (Springer
LNCS version isn't redistributable; author preprint is). Recommended
embedding the SVG and linking to Damien Stehlé's self-hosted preprint
instead of bundling the PDF.

**You**: Yes.
**Me**: Downloaded the SVG into `docs/src/assets/`, embedded it in the
tutorial and the algorithm explanation page (with attribution caption,
though public domain doesn't require it). Added Wikipedia + preprint
links to `explanation/algorithm.md` (References), `explanation/context.md`
(Further reading), and the `shortenW_in_UVW` docstring. Committed as
`8539901`, pushed.

---

## 9. v2.0.1 release — the stale-landing-page fix

**You**: The public docs site's landing page only links to the reference.
**Me**: Investigated. Found three layers of the problem, in order of
discovery:
1. Root `index.html` on `gh-pages` redirected to `stable/`, but `stable`
   was a symlink to `v1.0.4` (from March 2021).
2. No `v2.0.0/` directory existed on gh-pages at all — the v2.0.0 tag
   was pushed when `Documentor.yml` was still empty, so no tag-driven
   deploy ever happened for it.
3. **The real cause**: GitHub Pages was configured to serve from
   `gh-pages:/docs/`, not the root of gh-pages. Every Documenter deploy
   had been writing to `gh-pages:/v*.*.*/`, `gh-pages:/dev/`,
   `gh-pages:/stable/` — all invisible to the Pages server.

Proposed: cut v2.0.1 (carries the post-v2.0.0 docs/CI work), also fix
the cosmetic `repo=` warning in `docs/make.jl`.

**You**: Execute, fix the warning.
**Me**: Bumped Project.toml to 2.0.1, switched `repo=` to
`Documenter.Remotes.GitHub(...)`, committed, pushed, tagged v2.0.1,
pushed tag. CI built and deployed v2.0.1 to `gh-pages/v2.0.1/`.

The landing page still showed stale content. Diagnosed the GitHub
Pages serving-path misconfiguration via `gh api /pages`, fixed it with
`gh api --method PUT /pages -f source[path]=/`, triggered a Pages
rebuild, and verified the live landing page now shows the full quadrant
navigation.

---

## 10. `gh` CLI installation

**You**: How do I install gh so I don't have to do the release manually?
**Me**: `brew install gh && gh auth login`, interactive browser flow,
~30s. Listed the benefits: authenticated API (5000/hour vs 60),
`gh release create`, `gh run watch`, etc. — would have avoided the
rate-limit issues we hit during API polling.

**You**: gh is installed and authorized.
**Me**: Ran `gh release create v2.0.1 ...` — release created directly,
no manual step. Then `gh run list` to verify CI.

---

## 11. CI flake + iteration-cap bump

**You**: (Observed during CI watch) Runtests on main failed on
`RandUnimodMat3(10) * noisy_A`.
**Me**: Diagnosed: pre-existing flaky test. `RandUnimodMat3(10)` can
produce transforms that, times near-singular `A`, occasionally need
more than 15 iterations. The v2.0.1 tag-triggered run passed on the
same commit, confirming flakiness.

**You**: Seeding is brittle because Julia RNGs can change. Maybe just
raise the cap.
**Me**: Agreed, but flagged that "50" was a guess.

**You**: Does raising to 50 make sense? 15 was motivated.
**Me**: Ran an empirical study. Scaling of `DeviousMat(n)` gives ≈0.3
iterations/bit; extrapolating to Int64 worst case = ~20 iters. Float64
off-by-one `floor()` corrections add up to ~10 more. 50,000 randomised
stress trials produced max 23. Proposed **29** (= ~20 + ~10 headroom),
with a multi-line comment in the source that shows the derivation.

**You**: Let's go for 29.
**Me**: Changed the cap, updated the comment, `DeviousMat` docstring,
tutorial, `explanation/algorithm.md`, and `explanation/precision.md`.
Re-ran 20,000 trials of the previously-flaky pattern — max 19
iterations, zero cap hits. Committed `dc7666d`, pushed; both workflows
green.

**You**: Did it all pass?
**Me**: Yes — Runtests and Documentation both succeeded on `dc7666d`,
including the macOS-aarch64 runner that had flaked under the old cap.

---

## Final state at end of session

- **Version**: 2.0.1, released and deployed.
- **Docs**: Live at https://glwhart.github.io/MinkowskiReduction.jl with
  the full Diátaxis quadrant structure (tutorial, how-to, reference,
  explanation).
- **CI**: Runtests + Documentation workflows green on all three OSes;
  doctests exercised on every PR and push.
- **API change**: `minkReduce` returns the integer transform matrix `P`
  (the Spacey-enabling piece).
- **Iteration cap**: 29, motivated by theory + empirical stress tests.
- **Test count**: 67,406 (up from 64,377 at start of session).
- **Test reliability**: no flakes under the new cap, verified across
  20,000 randomised trials.

Outstanding items of note (not blockers, not worked on):

- The `.claude/` session-state directory appeared untracked in `git
  status`; adding it to `.gitignore` is a one-line follow-up.
- The v2.0.1 release has not been registered in the Julia General
  Registry (no `@JuliaRegistrator register` comment posted on the
  release).
