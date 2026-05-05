# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Run all tests
julia --project -e 'using Pkg; Pkg.test()'

# Run tests directly (faster iteration)
julia --project test/runtests.jl

# Build docs locally
cd docs && julia --project make.jl
```

## Architecture

This is a pure Julia package (no external dependencies beyond stdlib) implementing Minkowski lattice basis reduction for 2D and 3D lattices. The entire implementation lives in a single file: `src/MinkowskiReduction.jl`.

**Core algorithm (`mink_reduce`):** Iteratively sorts basis vectors by norm, then calls `shorten_w_in_uvw` to reduce the longest vector. That helper first calls `gauss_reduce` on the two shorter vectors (Euclidean-algorithm style 2D reduction), then projects the third vector onto the U-V parallelogram and subtracts integer multiples to minimize its distance to the origin. Repeats until convergence.

**Key exported functions:**
- `mink_reduce(U, V, W)` / `mink_reduce(M)` — main entry point; returns `(U, V, W, iterations)`
- `gauss_reduce(U, V)` — 2D Minkowski/Gauss reduction
- `is_mink_reduced(U, V, W)` / `is_mink_reduced(M)` — validates all 12 Minkowski conditions
- `orthogonality_defect(a, b, c)` / `orthogonality_defect(M)` — ratio of vector norms to `|det|`; equals 1 for orthogonal bases
- `rand_unimod_mat3`, `rand_unimod_mat2`, `devious_mat`, `is_permutation_matrix` — test utilities

**Both vector-tuple and matrix interfaces exist** for the main functions. The matrix form takes/returns columns as basis vectors.

**Test suite structure (`test/runtests.jl`):** Tests basic reduction, determinant preservation, `is_mink_reduced` (11 negative cases + positive cases), numerical robustness on FCC/BCC lattices with random noise across logarithmically-spaced noise levels, and edge cases (linearly dependent vectors throw `ErrorException`).

**`devious_mat(n)`** generates pathological unimodular matrices using Pisot numbers `(2+√3)^n` — used for stress-testing iteration counts.
