module MinkowskiReduction

using LinearAlgebra, Random
export GaussReduce, RandUnimodMat, RandLowerTri, minkReduce, DeviousMat, isMinkReduced, orthogonalityDefect
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
    # formed by the new U, V
    U, V = GaussReduce(U,V)
    # Project W into the U-V plane, call it T
    # n̂ is a unit vector ⟂ to U-V plane. Subtract multiples of this
    # from W to move it into the U-V plane. Dot product gives number of multiples. 
    n̂ = (U×V)/norm(U×V)
    T = W - W⋅n̂ * n̂
    # find multiples of U, V that move T inside the parallelogram formed by U, V 
    # (See notes from Rod Forcade in the support folder.)
    denom = (U⋅U)*(V⋅V)-(U⋅V)^2
    a = floor(((U⋅W)*(V⋅V)-(V⋅W)*(U⋅V))/denom)
    b = floor(((V⋅W)*(U⋅U)-(U⋅W)*(U⋅V))/denom)
    # Find the corner of the U-V cell closest to shifted T and shift T to be closest to origin
    W = W - a*U - b*V # Try the corner in "first quadrant"

    # Now try the other three corners, keep the shortest that is found
    norm(W) > norm(W - U) && (W = W - U)
    norm(W) > norm(W - V) && (W = W - V)
    norm(W) > norm(W - U - V) && (W = W - U - V)
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
    RandUnimodMat(n)

Generate a random unimodular 2x2 matrix. `n` is a small integer (number of row and column operations).

See also: `RandLowerTri(n)`, `FibonacciMat(n)`, `DeviousMat(n)`
"""
function RandUnimodMat(n)
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

end

