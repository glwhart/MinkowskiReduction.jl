module MinkowskiReduction

using LinearAlgebra, Random
export GaussReduce, RandUnimodMat2, RandLowerTri, minkReduce, DeviousMat, isMinkReduced, orthogonalityDefect, RandUnimodMat3, isPermutationMatrix
"""
    minkReduce(U, V, W)

Find the shortest equivalent basis of the lattice formed by {`U`, `V`, `W`}.

Returns a 4-tuple `(U′, V′, W′, n)`: the reduced basis vectors (with
`norm(U′) ≤ norm(V′) ≤ norm(W′)`) and the number of outer-loop iterations
`n` that were required. The reduction is not unique — see `research.md`
for discussion of sign, permutation, and boundary ambiguities.

```jldoctest
julia> U = [1, 2, 3]; V = [-1, 2, 3]; W = [3, 0, 4]; minkReduce(U,V,W)
([-2.0, 0.0, 0.0], [0.0, -2.0, 1.0], [-1.0, 2.0, 3.0], 2)
```
"""
function minkReduce(U, V, W)
    i = 0
    while true
        i +=1
        norms = [norm(U), norm(V), norm(W)]
        p = sortperm(norms)
        U,V,W = (U,V,W)[p] # sort into ascending order
        U,V,W = shortenW_in_UVW(U, V, W)
        i > 15 && error("minkReduce: Too many iterations")
#        println(U,V,W,i,"det",det(hcat(U,V,W)))
        norm(W) ≥ norm(V) ≥ norm(U) && break
    end
    #if !debug return U, V, W end
    return U, V, W, i
end

"""
    minkReduce(M)

Find the shortest equivalent basis of the lattice formed by the columns of
the 3×3 matrix `M`. Returns a 3×3 matrix whose columns are the reduced
basis vectors, in ascending norm order.

This is a convenience wrapper around [`minkReduce(U,V,W)`](@ref) that takes
and returns a matrix. The iteration count returned by the three-vector
form is discarded; call the three-vector form directly if you need it.

# Examples
```jldoctest
julia> minkReduce([1 -1 3; 2 2 0; 3 3 4])
3×3 Matrix{Float64}:
 -2.0   0.0  -1.0
  0.0  -2.0   2.0
  0.0   1.0   3.0
```
"""
function minkReduce(M)
    U,V,W = minkReduce(M[:,1],M[:,2],M[:,3])
    return hcat(U,V,W)
end

"""
    shortenW_in_UVW(U, V, W)

Shorten `W` by adding integer combinations of `U` and `V` to it, producing
the lattice vector of the coset `W + ℤU + ℤV` that is closest to the
origin. This is one inner step of the 3D Minkowski reduction algorithm of
Nguyen and Stehlé.

Returns the 3-tuple `(U′, V′, W′)`, where:
- `U′, V′` are the Gauss-reduced form of the input `(U, V)` pair (the
  function calls [`GaussReduce`](@ref) internally; the caller does not
  need to pre-reduce them),
- `W′ = W - aU′ - bV′` for integers `a, b` chosen so that `W′` is the
  shortest vector in that coset.

The geometric picture: after Gauss reduction, the (U′, V′)-parallelogram
contains the Voronoi cell of the 2D sublattice, so the closest lattice
point to the projection of `W` onto the U-V plane is one of the four
parallelogram corners. The function evaluates all four candidates and
picks the best.

Reference: Nguyen and Stehlé, *Low-Dimensional Lattice Basis Reduction
Revisited*, in Algorithmic Number Theory — ANTS VI (2004), LNCS 3076,
pp. 338–357.
"""
function shortenW_in_UVW(U,V,W)
    # If U, V are themselves mink reduced, then the projection of W (shifted by multiples of U,V)
    # that is the closest to the origin, will be contained in the parallelogram
    # formed by the new U, V. Find the corner of the parallelogram closest to W. Pick the vector from that corner to W as the new W.
    U, V = GaussReduce(U,V)
    # find multiples of U, V that move the projection of W inside the parallelogram formed by U, V
    denom = (U⋅U)*(V⋅V)-(U⋅V)^2
    a = floor(((U⋅W)*(V⋅V)-(V⋅W)*(U⋅V))/denom)
    b = floor(((V⋅W)*(U⋅U)-(U⋅W)*(U⋅V))/denom)
    W = W - a*U - b*V # Try the corner in "first quadrant"

    # Now find the corner of the parallelogram closest to the projection of W. Pick the vector from that corner to W as the new W.
    dists = [norm(W), norm(U - W), norm(V - W), norm(U + V - W)]
    case = argmin(dists)
    if case == 1 # Case 1 isn't necessary, but makes the logic clear
        W = W
    elseif case == 2
        W = W - U
    elseif case == 3
        W = W - V
    elseif case == 4
        W = W - U - V
    end
    return U, V, W
end

"""
    GaussReduce(U, V)

Reduce the 2D basis vectors {`U`, `V`} to the shortest equivalent basis,
using the classical Gauss–Lagrange algorithm (iterated Euclidean-style
reduction).

Returns a 2-tuple of reduced vectors in ascending norm order: if
`(a, b) = GaussReduce(U, V)` then `norm(a) ≤ norm(b)` (up to floating-
point tolerance). The pair satisfies the 2D Minkowski conditions
`norm(a) ≤ norm(b) ≤ norm(b ± a)`.

Throws an error if the input vectors are linearly dependent (parallel),
which is detected as NaN norms arising from the recursion's division by
zero.

# Examples
```jldoctest
julia> GaussReduce([5 8], [8 13])
([0.0 -1.0], [-1.0 0.0])
```
"""
function GaussReduce(U, V)
    maxval = max(abs.(U)...,abs.(V)...)
    if norm(U) > norm(V) U, V = V, U end
    i = 0
    while true
        V, U = U, V - round((U⋅V)/(U⋅U))*U
        i += 1
        if norm(U) > norm(V) || norm(U)≈norm(V) break; end
        # NaN norms arise when U and V are (numerically) parallel. This means the input
        # vectors are linearly *dependent*.
        isnan(norm(U)) && error("GaussReduce: input vectors are linearly dependent")
        i > 50 && error("GaussReduce: Too many iterations") # failsafe to break out if not converging
    end
    return V, U
end

"""
    orthogonalityDefect(a, b, c)

Compute the orthogonality defect of three basis vectors.

# Examples
```jldoctest
julia> orthogonalityDefect([1,1,0],[1,0,1],[0,1,1])
1.4142135623730954
```
"""
function orthogonalityDefect(a, b, c)
    return prod(norm.([a,b,c]))/abs((a×b)⋅c)
end

"""
    orthogonalityDefect(M)

Compute the orthogonality defect of the three column vectors of matrix `M`.

# Examples
```jldoctest
julia> M = [1 1 0; 1 0 1; 0 1 1]
3×3 Matrix{Int64}:
 1  1  0
 1  0  1
 0  1  1

julia> orthogonalityDefect(M)
1.4142135623730954
```
"""
function orthogonalityDefect(M::AbstractMatrix)
    size(M, 2) == 3 || error("Matrix must have exactly 3 columns")
    return orthogonalityDefect(M[:, 1], M[:, 2], M[:, 3])
end

"""
    RandUnimodMat2(n)

Generate a random unimodular (determinant `±1`) 2×2 integer matrix by
composing `n` random lower-triangular and upper-triangular integer shears.

Larger `n` produces matrices with larger entries and (in combination with
[`minkReduce`](@ref)) provides stress-test inputs for the 2D Gauss
reduction.

See also: [`RandLowerTri`](@ref), [`FibonacciMat`](@ref),
[`DeviousMat`](@ref), [`RandUnimodMat3`](@ref).
"""
function RandUnimodMat2(n)
    mat = RandLowerTri(1)
    for i ∈ 1:n
        mat = mat*RandLowerTri(1)
        mat = mat*transpose(RandLowerTri(1))
    end
    return mat
end

"""
    RandLowerTri(n)

Generate a 2×2 lower-triangular integer shear of the form

    [1 0;
     k 1]

where `k` is drawn uniformly at random from `-n:n`. The result is
unimodular (determinant = 1) by construction; it is used as a building
block for [`RandUnimodMat2`](@ref).

See also: [`RandUnimodMat2`](@ref), [`FibonacciMat`](@ref),
[`DeviousMat`](@ref).
"""
function RandLowerTri(n)
    return [1 0; rand(-n:n) 1]
end

"""
    FibonacciMat(k)

Generate a 2×2 matrix `[f2 f3; f1 f2]` whose entries are three
consecutive Fibonacci numbers (approximated via the closed-form Binet
expression with φ = 1.618…). The resulting matrix is unimodular and
deliberately ill-conditioned: its columns are nearly parallel, which
makes it a hard stress-test for 2D Gauss reduction.

`k` selects which Fibonacci triple is used; larger `k` gives larger
entries and more extreme ill-conditioning. The function errors on
`Int64` overflow, which occurs around `k ≈ 92`.

See also: [`RandUnimodMat2`](@ref), [`DeviousMat`](@ref).
"""
function FibonacciMat(k)
    f1 = round(Int64,1.61803398875^k/sqrt(5))
    f2 = round(Int64,1.61803398875^(k+1)/sqrt(5))
    f3 = f1 + f2
    any(i -> i < 1, [f1 f2 f3]) && error("Overflow in FibonacciMat function")
    return [f2 f3; f1 f2]
end

"""
    DeviousMat(n)

Generate a unimodular 3×3 integer matrix that is a heavily disguised
simple-cubic basis and requires a large number of reduction steps to
recognize as such. It is the hardest known 3D stress test for the
Nguyen–Stehlé greedy algorithm.

The matrix entries grow like the Pisot number `(2+√3)ⁿ`, so the cost of
reducing it scales linearly with `n`. Because the entries are stored as
`Int64`, `n` must satisfy `3 ≤ n ≤ 26`; `n = 27` silently overflows.
`n = 26` is the largest representable instance and requires exactly 15
outer iterations of [`minkReduce`](@ref) — the value used to set the
iteration cap in `minkReduce`.

(Construction due to Rod Forcade, private communication, Feb 1 2020.)

See also: [`RandUnimodMat3`](@ref), [`FibonacciMat`](@ref).
"""
function DeviousMat(n)
    n < 3 && error("for DeviousMat, n > 2")
    u,v = round(Int64,(2+√3)^n/(2*√3)), round(Int64,(2+√3)^n/2)
    a,b = convert(Int64,(u+v+1)/2), -u
    c,d = a-1, v-u
    return [a b c; b d b; c b a]
end

"""
    isMinkReduced(U,V,W)

Check if the basis {`U`,`V`,`W`} is Minkowski reduced.

Each of the 12 defining inequalities is checked up to a floating-point
tolerance that scales with the largest of ‖U‖, ‖V‖, or ‖W‖.
"""
function isMinkReduced(U,V,W)
    tol = eps(max(norm(U), norm(V), norm(W)))
    if norm(U) > norm(V)+tol println("Condition 1 failed"); return false end
    if norm(V) > norm(W)+tol println("Condition 2 failed"); return false end
    if norm(V) > norm(U+V)+tol println("Condition 3 failed"); return false end
    if norm(V) > norm(U-V)+tol println("Condition 4 failed"); return false end
    if norm(W) > norm(U+W)+tol println("Condition 5 failed"); return false end
    if norm(W) > norm(U-W)+tol println("Condition 6 failed"); return false end
    if norm(W) > norm(V+W)+tol println("Condition 7 failed"); return false end
    if norm(W) > norm(V-W)+tol println("Condition 8 failed"); return false end
    if norm(W) > norm(U+V+W)+tol println("Condition 9 failed"); return false end
    if norm(W) > norm(U-V+W)+tol println("Condition 10 failed"); return false end
    if norm(W) > norm(U+V-W)+tol println("Condition 11 failed"); return false end
    if norm(W) > norm(U-V-W)+tol println("Condition 12 failed"); return false end
    return true
end

"""
    isMinkReduced(M)

Check if the basis formed by the columns of the 3×3 matrix `M` is
Minkowski reduced. Convenience wrapper around
[`isMinkReduced(U,V,W)`](@ref).
"""
function isMinkReduced(M)
    return isMinkReduced(M[:,1],M[:,2],M[:,3])
end

"""
    RandUnimodMat3(k=10)

Generate a random `3×3` unimodular matrix (determinant `±1`).

The algorithm starts with the identity matrix and performs `k` random
integer elementary operations that preserve the determinant:

1. Add an integer multiple of one *row* to another row.
2. Optionally, add an integer multiple of one *column* to another column.

Because these elementary operations have determinant `1`, the resulting
matrix is guaranteed to be unimodular.  The optional column operations are
included only to improve the "randomness" of the output.

`k` controls how many elementary operations are applied (default `10`).
Larger values of `k` will typically lead to matrices with larger
(integer-sized) entries.

# Examples
```jldoctest
julia> M = RandUnimodMat3()
3×3 Matrix{Int64}:
  1   0  -3
 -5   1  29
  4   0 -23

julia> det(M)
1
```
"""
function RandUnimodMat3(k::Integer = 10)
    k < 1 && error("RandUnimodMat3: k must be positive")
    M = Matrix{Int64}(I, 3, 3)

    for _ in 1:k
        # ----- Row shear -----
        r1, r2 = rand(1:3, 2)
        while r2 == r1
            r2 = rand(1:3)
        end
        m = rand(-k:k)
        m == 0 && (m = 1) # ensure a non-zero multiple
        @inbounds M[r1, :] .+= m .* M[r2, :]

        # ----- (Optional) Column shear -----
        if rand(Bool)
            c1, c2 = rand(1:3, 2)
            while c2 == c1
                c2 = rand(1:3)
            end
            n = rand(-k:k)
            n == 0 && (n = 1)
            @inbounds M[:, c1] .+= n .* M[:, c2]
        end
    end

    # After the shear operations det(M) is exactly 1.  Randomly flip the sign
    # of one row to obtain determinant ±1 with equal probability.
    if rand(Bool)
        row = rand(1:3)
        M[row, :] .*= -1
    end
    return M
end

"""
    isPermutationMatrix(M; atol=√eps()) -> Bool

Return `true` if the 3×3 matrix `M` is a (signed) permutation of the identity,
i.e. each row and each column contains exactly one entry whose absolute value
is ≈ 1 and all other entries are ≈ 0.  The optional keyword `atol` sets the
absolute tolerance used by `isapprox` when comparing entries to 0/1.

A *signed* permutation matrix differs from a permutation matrix only by having
rows possibly multiplied by –1, so `det(M) = ±1`.

# Examples
```jldoctest
julia> isPermutationMatrix(I)
true

julia> isPermutationMatrix([0 1 0; 0 0 1; 1 0 0])  # a permutation matrix
true

julia> isPermutationMatrix([0 -1 0; 0 0 1; -1 0 0]) # signed permutation
true
```
"""
function isPermutationMatrix(M::AbstractMatrix{<:Real}; atol = sqrt(eps()))
    # Must be 3×3
    size(M) == (3,3) || return false

    # Check structure: exactly one (≈) ±1 per row/column, rest ≈0
    for i in 1:3
        rowvals = M[i, :]
        colvals = M[:, i]
        # Identify entries whose |val| ≈ 1
        ones_in_row = count(i -> isapprox(abs(rowvals[i]), 1; atol=atol), 1:3)
        zeros_in_row = count(i -> isapprox(rowvals[i], 0; atol=atol), 1:3)
        ones_in_col = count(i -> isapprox(abs(colvals[i]), 1; atol=atol), 1:3)
        zeros_in_col = count(i -> isapprox(colvals[i], 0; atol=atol), 1:3)
        if !(ones_in_row == 1 && zeros_in_row == 2 && ones_in_col == 1 && zeros_in_col == 2)
            return false
        end
    end
    return true
end

end
