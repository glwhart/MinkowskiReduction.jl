# MinkowskiReduction.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://glwhart.github.io/MinkowskiReduction.jl)
[![Build Status](https://travis-ci.com/glwhart/MinkowskiReduction.jl.svg?branch=main)](https://travis-ci.com/glwhart/MinkowskiReduction.jl)
[![Coverage](https://coveralls.io/repos/github/glwhart/MinkowskiReduction.jl/badge.svg?branch=main)](https://coveralls.io/github/glwhart/MinkowskiReduction.jl?branch=main)

Reduces a basis for a three-dimensional lattice to the basis with the shortest possible basis vectors. Equivalently, the code finds the basis that is as close to orthogonal as possible. This is known as Minkowski reduction. (See _Geometrie der Zahlen_ Hermann Minkowski 1910). In higher dimensions, the famous [LLL lattice reduction](https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm) is commonly used. Lattice reduction is NP-complete in d > 4 dimensions, but a polynomial time algorithm exists for three and four dimensions.[Phong Nguyen and Damien Stehlé, "Low-Dimensional Lattice Basis Reduction Revisited
"](https://link.springer.com/chapter/10.1007/978-3-540-24847-7_26) This implementation for three dimensions was inspired by this work.

This code is useful in Density Functional Theory calculations to convert a crystal structure to a compact basis best suited for accurate calculations. It is also useful for reducing a lattice prior to calculating the pointgroup because the symmetry of a reduced lattice can be found in a fixed (relatively small) number of steps. (See [Pointgroup finder Pointy.jl](https://github.com/glwhart/Pointy.jl)

# Example 1: Reduce a horribly skew basis
```
julia> bigM = DeviousMat(26) # Matrix whose columns are an extremely skew basis for the lattice of integers
3×3 Array{Int64,2}:
  292755045568446  -214311567528244   292755045568445
 -214311567528244   156886956080403  -214311567528244
  292755045568445  -214311567528244   292755045568446
  
julia> U,V,W = mapslices(x->[x], bigM, dims=2) # Grab columns for input to Minkowski reduction routine
3×1 Array{Array{Int64,1},2}:
 [292755045568446, -214311567528244, 292755045568445]
 [-214311567528244, 156886956080403, -214311567528244]
 [292755045568445, -214311567528244, 292755045568446]
  
  julia> minkReduce(U,V,W)
  ([-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0])
  ```
# Example 2: See how many steps are required to reduce a basis
```
julia> U = [292755045568446, -214311567528244, 292755045568445]; 
julia> V = [-214311567528244, 156886956080403, -214311567528244];
julia> W = [292755045568445, -214311567528244, 292755045568446]

# Use an optional argument to trigger the 4th item in the returned tuple
# to indicate the number of steps required to reduce the lattice
julia> minkReduce(U,V,W,true)
([-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0], 15)
```

