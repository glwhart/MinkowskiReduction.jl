# MinkowskiReduction.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://glwhart.github.io/MinkowskiReduction.jl)
[![Runtests](https://github.com/glwhart/MinkowskiReduction.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/glwhart/MinkowskiReduction.jl/actions/workflows/Runtests.yml)
[![codecov](https://codecov.io/gh/glwhart/MinkowskiReduction.jl/branch/main/graph/badge.svg?token=2SG9ZXKC2W)](https://codecov.io/gh/glwhart/MinkowskiReduction.jl)

Reduces a basis for a three-dimensional lattice to the basis with the shortest possible basis vectors. Equivalently, the code finds the basis that is as close to orthogonal as possible. This is known as Minkowski reduction. (See _Geometrie der Zahlen_ Hermann Minkowski 1910). In higher dimensions, the famous [LLL lattice reduction](https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm) is commonly used. Lattice reduction is NP-complete in d > 4 dimensions, but a polynomial time algorithm exists for three and four dimensions. (See [Phong Nguyen and Damien Stehlé, "Low-Dimensional Lattice Basis Reduction Revisited
"](https://link.springer.com/chapter/10.1007/978-3-540-24847-7_26)) The implementation in this repo for three dimensions (and two dimensions) was inspired by this work.

This code is useful in Density Functional Theory calculations to convert a crystal structure to a compact basis best suited for accurate calculations. It is also useful for reducing a lattice prior to calculating the pointgroup because the symmetry of a reduced lattice can be found in a fixed (relatively small) number of steps. (See pointgroup and spacegroup calculator [Spacey.jl](https://github.com/glwhart/Spacey.jl))


# Example 1: Reduce a slightly skew basis for a simple cubic basis
Consider a simple cubic lattice. The "natural" basis is just the _standard basis_: (1,0,0), (0,1,0), and (0,0,1) and any equivalent (but skew) basis should reduce to this one. As an example, take this basis: (1,0,0), (0,1,0), and (1,1,1). The first two are orthogonal but the third one is not. Reduce it with the `minkReduce` function.
```
julia> a1 = [1.,0.,0.]
julia> a2 = [0.,1.,0.]
julia> a3 = [1.,1.,1.]
julia> minkReduce(a1,a2,a3)
([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0])
```
The output is the standard basis, as expected.

# Example 2: Reduce a skew basis for a hexagonal lattice
```
julia> a1 = [1.,0.,0.]
julia> a2 = [1.5,√3/2,0.]
julia> a3 = [0,0.,√(8/3)]
julia> minkReduce(a1,a2,a3)
([0.5, 0.8660254037844386, 0.0], [0.5, -0.8660254037844386, 0.0], [0.0, 0.0, 1.632993161855452])
```
# Example 3: Reduce a horribly skew basis
```
julia> bigM = DeviousMat(26) # Matrix whose columns are an extremely skew basis for a simple cubic lattice
3×3 Array{Int64,2}:
  292755045568446  -214311567528244   292755045568445
 -214311567528244   156886956080403  -214311567528244
  292755045568445  -214311567528244   292755045568446
  
julia> U,V,W = mapslices(x->[x], bigM, dims=2) # Grab columns for input to Minkowski reduction routine
3×1 Array{Array{Int64,1},2}:
 [292755045568446, -214311567528244, 292755045568445]
 [-214311567528244, 156886956080403, -214311567528244]
 [292755045568445, -214311567528244, 292755045568446]
  
  julia> u,v,w = minkReduce(U,V,W) # Reduce the basis and store in new variables
  ([-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0], 15)
  ```
  The output is the standard basis again (modulo negative signs and ordering). The fourth return value (an integer) is the number of steps required to reduce the basis. In this extreme case, 15 steps are required.
