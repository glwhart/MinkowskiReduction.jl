module MinkowskiReduction

using LinearAlgebra, Random
export GaussReduce, RandUnimodMat2, RandLowerTri, minkReduce, DeviousMat, isMinkReduced, orthogonalityDefect, RandUnimodMat3, isPermutationMatrix
"""
    minkReduce(U, V, W, debug=false)

Find the shortest equivalent basis of that lattice formed by {`U`, `V`, `W`}

```jldoctest
julia> U = [1, 2, 3]; V = [-1, 2, 3]; W = [3, 0, 4]; minkReduce(U,V,W)
([-2.0, 0.0, 0.0], [0.0, -2.0, 1.0], [-1.0, 2.0, 3.0])
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

""" minkReduce(M) 

Find the shortest equivalent basis of that lattice formed by the columns of `M`."""
function minkReduce(M)
    U,V,W = minkReduce(M[:,1],M[:,2],M[:,3])
    return hcat(U,V,W)
end

"""
    shortenW_in_UVW

Reduce vector W so that it is as close to the origin as possible.

Subtract multiples of U and V from W. W will remain in an affine plane, which 
is parallel to the U-V plane but which passes through the end of the W vector. 
(See Lecture notes in computer science, ISSN 0302-974, ANTS - VI : algorithmic 
number theory, 2004, vol. 3076, pp. 338-357 ISBN 3-540-22156-5)
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

Reduce the basis vectors {`U`, `V`} to the shortest possible basis.

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
        isnan(norm(U)) && error("GaussReduce: input vectors were linearly independent") 
        i > 50 && error("GaussReduce: Too many iterations") # failsafe to break out if not converging
    end
    return V, U
end

"""
    orthogonalityDefect(a,b,c)

Compute the orthogonality defect of three basis vectors.

# Examples
```jldoctest
julia> orthogonalityDefect([1,1,0],[1,0,1],[0,1,1])
1.4142135623730954
```
"""
function orthogonalityDefect(a,b,c)
    return prod(norm.([a,b,c]))/abs((a×b)⋅c)
end

"""
    RandUnimodMat2(n)

Generate a random unimodular 2x2 matrix. `n` is a small integer (number of row and column operations).

See also: `RandLowerTri(n)`, `FibonacciMat(n)`, `DeviousMat(n)`
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
Generate a random 2x2 matrix of the form [1 0; 0 ±n].

See also: `RandUnimodMat(n)`, `FibonacciMat(n)`, `DeviousMat(n)`
"""
function RandLowerTri(n)
    return [1 0; rand(-n:n) 1]
end

"""
    FibonacciMat(k)
Generate a 2x2 matrix of the form [f2 f3; f1 f2] where f1, f2, f3 are consecutive Fibonacci-like numbers

See also: `RandUnimodMat(n)` and `DeviousMat(n)`

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

Make a unimodular 3x3 matrix that requires a large number of steps to reduce
(See email from Rod Feb 1 2020)

`n` dictates the size of the entries
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

Check if the basis {`U`,`V`,`W`} is Minkoswki reduced.
    
"""
function isMinkReduced(U,V,W)
    if norm(U) > norm(V)+eps() println("Condition 1 failed"); return false end
    if norm(V) > norm(W)+eps() println("Condition 2 failed"); return false end
    if norm(V) > norm(U+V)+eps() println("Condition 3 failed"); return false end
    if norm(V) > norm(U-V)+eps() println("Condition 4 failed"); return false end
    if norm(W) > norm(U+W)+eps() println("Condition 5 failed"); return false end
    if norm(W) > norm(U-W)+eps() println("Condition 6 failed"); return false end
    if norm(W) > norm(V+W)+eps() println("Condition 7 failed"); return false end
    if norm(W) > norm(V-W)+eps() println("Condition 8 failed"); return false end
    if norm(W) > norm(U+V+W)+eps() println("Condition 9 failed"); return false end
    if norm(W) > norm(U-V+W)+eps() println("Condition 10 failed"); return false end
    if norm(W) > norm(U+V-W)+eps() println("Condition 11 failed"); return false end
    if norm(W) > norm(U-V-W)+eps() println("Condition 12 failed"); return false end
    return true
end

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

