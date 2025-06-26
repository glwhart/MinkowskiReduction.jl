using MinkowskiReduction
using Test
using Random
using LinearAlgebra

@testset "MinkowskiReduction.jl" begin
    U=[1, 2, 3];V=[-1, 2, 3];W=[3, 0, 4]
    @test all(minkReduce(U,V,W) .≈ ([-2.0, 0.0, 0.0], [0.0, -2.0, 1.0], [-1.0, 2.0, 3.0],2))
    @test orthogonalityDefect(minkReduce(U,V,W)[1:3]...)≈1.0458250331675945
    @test orthogonalityDefect(U,V,W)==4.375
    fcc = [[1,1,0],[1,0,1],[0,1,1]]
    @test orthogonalityDefect(fcc...)≈1.4142135623730954
    m = DeviousMat(26) # Largest size that doesn't overflow
    @test abs(det(hcat(minkReduce(m[:,1],m[:,2],m[:,3])[1:3]...)))==1
    @test all(MinkowskiReduction.shortenW_in_UVW(U,V,W) .≈ ([-2.0, 0.0, 0.0], [1.0, 2.0, 3.0], [0.0, -2.0, 1.0]))
    @test (m = DeviousMat(26); all(minkReduce(m[:,1],m[:,2],m[:,3]) .≈ ([1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0], 15)))                                                                  
    @test (m = DeviousMat(20); all(minkReduce(m[:,1],m[:,2],m[:,3]) .≈ ([0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 12)))
    # These next few are not robust, but testing something that depends on the random number generator is difficult
    @test det(MinkowskiReduction.RandLowerTri(35))≈1.0 
    @test det(MinkowskiReduction.RandUnimodMat2(4))≈1.0 
    @test det(MinkowskiReduction.RandLowerTri(45))≈1.0 
    @test det(MinkowskiReduction.RandUnimodMat2(6))≈1.0 
    @test (MinkowskiReduction.FibonacciMat(45)==[1836311903 2971215073; 1134903170 1836311903])
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
    U = -[2,0,0]; V = [0,2,0]; W=U+V # Linearly dependent basis
    @test_throws ErrorException minkReduce(U,V,W)
    @test_throws ErrorException minkReduce(hcat(U,V,W))
    U = -[2.,0,0]; V = [0,2.,0]; W=U+V+[0,0,1e-170]
    @test_throws ErrorException minkReduce(U,V,W)
    @test_throws ErrorException minkReduce(hcat(U,V,W))
    @test (U = -[2.,0,0]; V = [0,2.,0]; W=U+V+[0,0,1e-150];minkReduce(U,V,W)==
    ([0.0, 0.0, 1.0e-150], [-2.0, -0.0, -0.0], [0.0, 2.0, 0.0],2))
    @test (U = -[2.,0,0]; V = [0,2.,0]; W=U+V+[0,0,1e-150];minkReduce(hcat(U,V,W))==
    hcat([0.0, 0.0, 1.0e-150], [-2.0, -0.0, -0.0], [0.0, 2.0, 0.0]))
    #U = -[2^63-1,0,0]; V = [0,2^63-1,0]; W=U+V+[0,0,1e-150]; 
    #@test_throws InexactError minkReduce(U,V,W) # silently overflows with later versions of julia
    for i ∈ 1:100
        @test isPermutationMatrix(vcat(shuffle!([[1 0 0],[0 1 0],[0 0 1]])...))
    end
    A = [-1.7233692904465637e-9 0.0025000286163305314 0.0024999524070139123; 0.0024999730832105326 2.591951474986415e-8 0.0025000013059337757; 0.0024999629647447772 0.002499967442175358 -1.891891458670977e-8]
    for a ∈ logrange(1e-15,1e15,30)
        @test isMinkReduced(minkReduce(A*a))
    end
    for i ∈ 1:10
        M = RandUnimodMat3(10)
        @test isMinkReduced(minkReduce(M*A))
    end
    A = [0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0]
    for ε ∈ logrange(1e-5,1e-1,10)
        for i ∈ 1:100
            noise = (2*rand(3,3).-1)*ε
            @test isMinkReduced(minkReduce(A+noise))
            @test isMinkReduced(minkReduce(hcat(eachcol(A+noise)...)))
        end
    end
    # Test of noise levesls for BCC
    A = [1.0 1.0 -1.0; 1.0 1.0 1.0; -1.0 1.0 1.0]
    for ε ∈ logrange(1e-5,1.1e-1,20)
        for i ∈ 1:100
            noise = (2*rand(3,3).-1)*ε
            @test isMinkReduced(minkReduce(A+noise))
        end
    end
    # Loop over aspect ratios and noise levels for BCC
    # With large aspect ratios, the code is not as robust to noise (but these are huge aspect ratios)
    for ε ∈ logrange(1e-5,1.1e-1,20)
        for ar ∈ logrange(1e-8,1.1e8,30)
            aspect_ratio = [1 0 0; 0 1 0; 0 0 ar]
            A = aspect_ratio*[1.0 1.0 -1.0; 1.0 1.0 1.0; -1.0 1.0 1.0]
            for i ∈ 1:50
                noise = (2*rand(3,3).-1)*ε
                @test isMinkReduced(minkReduce(A+noise))
            end
        end
    end
    # Slightly more sensitive to noise for FCC (might just be due to different lattice parameter)
    for ε ∈ logrange(1e-5,1.1e-1,20)
        for ar ∈ logrange(1e-7,1.1e7,30)
            aspect_ratio = [1 0 0; 0 1 0; 0 0 ar]
            A = aspect_ratio*[0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0]
            for i ∈ 1:50
                noise = (2*rand(3,3).-1)*ε
                @test isMinkReduced(minkReduce(A+noise))
            end
        end
    end
end




