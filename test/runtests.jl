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
    #@test (Random.seed!(1234); MinkowskiReduction.RandLowerTri(35)==[1 0;-28 1])
    @test (MinkowskiReduction.FibonacciMat(45)==[1836311903 2971215073; 1134903170 1836311903])
    #@test (Random.seed!(1234);MinkowskiReduction.RandUnimodMat(4)==[-1 -1; 0 -1])
    #@test (Random.seed!(1234);MinkowskiReduction.RandUnimodMat(2)==[0 1; -1 -1])
    @test !(U = [1, 2, 9]; V = [-3, 4, 3]; W = [3, 0, 4]; isMinkReduced(U,V,W))
    @test !(U = [1, 2, 3]; V = [-3, 4, 3]; W = [3, 0, 4]; isMinkReduced(U,V,W))
    @test !(U = [1, 0, 0]; V = [-1, 0, 1]; W = [3, 0, 4]; isMinkReduced(U,V,W))
    @test !(U = [1, 2, 3]; V = [-1, 2, 3]; W = [3, 0, 4]; isMinkReduced(U,V,W))
    @test !(U = [1, 0, 0]; V = [0,-4,-3]; W = [-3, 0, 4]; isMinkReduced(U,V,W))
    @test !(U = [1, 0, 0]; V = [0, 4, 3]; W = [3, 0, 4]; isMinkReduced(U,V,W))
    @test !(U = [0, 4, 3]; V = [-5, 0, 0]; W = [3, 0, 4]; isMinkReduced(U,V,W))
    @test !(U = [1,0,0]; V = [0,1,0]; W = [0,1,-2]; isMinkReduced(U,V,W))
    @test !(U = [1,0,-2]; V = [-1,1,-2]; W = [0,0,3]; isMinkReduced(U,V,W))
    @test !(δ=.50;U = [1,0,0]; V = [-.4,1,0]; W = [-δ,-δ,5]; isMinkReduced(U,V,W))
    @test !(U = [2,0,0]; V = [-.9,2,0]; W = [-1,-1,10]; isMinkReduced(U,V,W))
    @test !(U = [2,0,0]; V = -[-.9,2,0]; W = [-1,-1,10]; isMinkReduced(U,V,W))
    @test !(U = -[2,0,0]; V = -[-.9,2,0]; W = [-1,-1,10]; isMinkReduced(U,V,W))
    @test !(U = -[2,0,0]; V = [-.9,2,0]; W = [-1,-1,10]; isMinkReduced(U,V,W))
    @test (U = -[2,0,0]; V = [0,2,0]; W = [-1,-1,10]; isMinkReduced(U,V,W))
    U = -[2,0,0]; V = [0,2,0]; W=U+V # Linearly dependedent basis
    @test_throws ErrorException minkReduce(U,V,W)
    U = -[2.,0,0]; V = [0,2.,0]; W=U+V+[0,0,1e-170]
    @test_throws ErrorException minkReduce(U,V,W)
    @test (U = -[2.,0,0]; V = [0,2.,0]; W=U+V+[0,0,1e-150];minkReduce(U,V,W)==
    ([0.0, 0.0, 1.0e-150], [-2.0, -0.0, -0.0], [0.0, 2.0, 0.0]))
    U = -[2^63-1,0,0]; V = [0,2^63-1,0]; W=U+V+[0,0,1e-150]; 
    @test_throws InexactError minkReduce(U,V,W)
end
