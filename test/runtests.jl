using MinkowskiReduction
using Test

@testset "MinkowskiReduction.jl" begin
    U=[1, 2, 3];V=[-1, 2, 3];W=[3, 0, 4]
    @test all(minkReduce(U,V,W) .≈ ([-2.0, 0.0, 0.0], [0.0, -2.0, 1.0], [-1.0, 2.0, 3.0]))
    @test orthogonalityDefect(minkReduce(U,V,W)...)≈1.0458250331675945
    @test orthogonalityDefect(U,V,W)==4.375
    fcc = [[1,1,0],[1,0,1],[0,1,1]]
    @test orthogonalityDefect(fcc...)≈1.4142135623730954
    m = DeviousMat(26) # Largest size that doesn't overflow
    @test det(hcat(minkReduce(m[:,1],m[:,2],m[:,3])...))==1
end
