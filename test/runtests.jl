using MinkowskiReduction
using Test
using Random
using LinearAlgebra

@testset "MinkowskiReduction.jl" begin
    U=[1, 2, 3];V=[-1, 2, 3];W=[3, 0, 4]
    @test all(minkReduce(U,V,W) .≈ ([-2.0, 0.0, 0.0], [0.0, -2.0, 1.0], [-1.0, 2.0, 3.0]))
    @test orthogonalityDefect(minkReduce(U,V,W)...)≈1.0458250331675945
    @test orthogonalityDefect(U,V,W)==4.375
    fcc = [[1,1,0],[1,0,1],[0,1,1]]
    @test orthogonalityDefect(fcc...)≈1.4142135623730954
    m = DeviousMat(26) # Largest size that doesn't overflow
    @test det(hcat(minkReduce(m[:,1],m[:,2],m[:,3])...))==1
    @test all(MinkowskiReduction.shortenW_in_UVW(U,V,W) .≈ ([-2.0, 0.0, 0.0], [1.0, 2.0, 3.0], [-2.0, -2.0, 1.0]))
    @test (m = DeviousMat(26); all(minkReduce(m[:,1],m[:,2],m[:,3],true) .≈ ([-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0], 15)))
    @test (m = DeviousMat(20); all(minkReduce(m[:,1],m[:,2],m[:,3],true) .≈ ([0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 12)))
    @test (Random.seed!(1234); MinkowskiReduction.RandLowerTri(35)==[1 0;-28 1])
    @test (MinkowskiReduction.FibonacciMat(45)==[1836311903 2971215073; 1134903170 1836311903])
end
