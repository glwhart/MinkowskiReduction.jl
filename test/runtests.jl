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
    @test orthogonalityDefect(hcat(fcc...))≈1.4142135623730954
    @test orthogonalityDefect(minkReduce(U,V,W)[1:3]...)≈1.0458250331675945
    @test orthogonalityDefect(hcat([U,V,W]...))==4.375
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
    for i ∈ 1:100
        @test !isPermutationMatrix(vcat(shuffle!([[1 0 1],[0 1 0],[0 0 1]])...))
        @test !isPermutationMatrix(vcat(shuffle!([[1 0 0],[0 0 0],[0 0 1]])...))
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

    # ========================================================================
    # Property-based tests. These do not verify specific outputs; they assert
    # invariants that any correct Minkowski reducer must satisfy, so they
    # catch regressions from future refactoring without being brittle against
    # the non-uniqueness of reduction.
    # ========================================================================

    seed_bases = [
        Matrix{Float64}(I, 3, 3),                        # simple cubic
        [0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0],         # FCC
        [-0.5 0.5 0.5; 0.5 -0.5 0.5; 0.5 0.5 -0.5],      # BCC
    ]

    # --- 1. minkReduce invariants on random unimodular transforms ------------
    # For each canonical Bravais lattice, apply random unimodular transforms
    # and check: the reduced output satisfies all 12 Minkowski conditions,
    # |det| is preserved, and the orthogonality defect does not increase.
    for A ∈ seed_bases
        for i ∈ 1:50
            # k=4 keeps entries of RandUnimodMat3 small enough that `det` of
            # the resulting matrix is accurate to machine precision; larger
            # k would exercise the reducer equally well but introduce
            # spurious cancellation error in the reference det itself.
            M = A * RandUnimodMat3(4)
            R = minkReduce(M)
            @test isMinkReduced(R)
            @test isapprox(abs(det(R)), abs(det(M)); rtol = 1e-10)
            @test orthogonalityDefect(R) ≤ orthogonalityDefect(M) * (1 + 1e-10)
        end
    end

    # --- 2. Idempotence ------------------------------------------------------
    # Reducing an already-reduced basis must be a no-op.
    for i ∈ 1:50
        B = minkReduce(Float64.(RandUnimodMat3(8)))
        @test minkReduce(B) ≈ B
    end

    # --- 3. GaussReduce (currently has no dedicated tests) -------------------
    # 2D Minkowski conditions on output, |det| preservation, idempotence,
    # and error on linearly dependent input.
    for i ∈ 1:100
        U = rand(-100:100, 2)
        V = rand(-100:100, 2)
        (iszero(U) || iszero(V)) && continue
        d = abs(det(hcat(U, V)))
        d == 0 && continue                   # skip linearly dependent inputs
        short, long = GaussReduce(U, V)      # GaussReduce returns (shorter, longer)
        @test norm(short) ≤ norm(long) + sqrt(eps())
        @test norm(long) ≤ norm(long + short) + sqrt(eps())
        @test norm(long) ≤ norm(long - short) + sqrt(eps())
        @test abs(det(hcat(short, long))) ≈ d
    end
    # Idempotence on the docstring example
    let (a1, b1) = GaussReduce([5, 8], [8, 13])
        a2, b2 = GaussReduce(a1, b1)
        @test a2 ≈ a1 && b2 ≈ b1
    end
    # Parallel vectors are detected and raise an error
    @test_throws ErrorException GaussReduce([1.0, 0, 0], [2.0, 0, 0])

    # --- 4. Sign and permutation invariance of the norm multiset -------------
    # Signed permutations of the input columns preserve the lattice, and so
    # must leave the sorted sequence of reduced norms unchanged (even though
    # the output column order / signs may differ due to non-uniqueness).
    norm_multiset(M) = sort(norm.(eachcol(M)))
    perms_3 = [[1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]]
    for A ∈ seed_bases
        ref = norm_multiset(minkReduce(A))
        for p ∈ perms_3, s1 ∈ (-1,1), s2 ∈ (-1,1), s3 ∈ (-1,1)
            Ap = A[:, p] .* reshape([s1, s2, s3], 1, 3)
            @test norm_multiset(minkReduce(Ap)) ≈ ref
        end
    end

    # --- 5. Scale invariance of the iteration count --------------------------
    # The algorithm makes only ratio-based decisions, so uniformly scaling
    # the basis must not change how many outer iterations it takes.
    U0, V0, W0 = [1.0, 2, 3], [-1.0, 2, 3], [3.0, 0, 4]
    _, _, _, n_ref = minkReduce(U0, V0, W0)
    for α ∈ logrange(1e-8, 1e8, 20)
        _, _, _, n = minkReduce(α*U0, α*V0, α*W0)
        @test n == n_ref
    end

    # --- 6. Canonical Bravais lattices reduce to known invariants ------------
    # For primitive cells of several Bravais lattice types with known
    # successive minima and orthogonality defect, apply random unimodular
    # transforms and verify that the reducer recovers those invariants.
    # We check sorted norms and the defect rather than specific basis
    # matrices, because Minkowski reduction is non-unique (sign,
    # permutation, and boundary-equality ambiguities).
    #
    # Rhombohedral primitive cell with vector length `a` and equal pairwise
    # angle α ∈ (0, 2π/3).
    rhombohedral_basis(a, α) = begin
        c, s = cos(α), sin(α)
        z = sqrt((1 + 2c) * (1 - c) / (1 + c))
        a .* [1.0  c            c;
              0.0  s            c*(1-c)/s;
              0.0  0.0          z]
    end
    # Orthogonality defect of the rhombohedral cell, derived from
    # det(Gram) = 1 − 3cos²α + 2cos³α for unit-length vectors.
    rhombohedral_defect(α) = 1 / sqrt(1 - 3*cos(α)^2 + 2*cos(α)^3)

    # Each entry: (name, primitive basis, expected sorted norms, expected defect)
    bravais = [
        ("Cubic-P",
            Matrix{Float64}(I, 3, 3),
            [1.0, 1.0, 1.0],     1.0),
        ("FCC",
            [0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0],
            fill(1/√2, 3),       √2),
        ("BCC",
            [-0.5 0.5 0.5; 0.5 -0.5 0.5; 0.5 0.5 -0.5],
            fill(√3/2, 3),       3√3/4),
        ("Tetragonal-P (c/a=1.5)",
            [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.5],
            [1.0, 1.0, 1.5],     1.0),
        ("Orthorhombic-P",
            [1.0 0.0 0.0; 0.0 1.3 0.0; 0.0 0.0 1.7],
            [1.0, 1.3, 1.7],     1.0),
        ("Hexagonal-P (c/a=1.6)",
            [1.0 -0.5 0.0; 0.0 √3/2 0.0; 0.0 0.0 1.6],
            [1.0, 1.0, 1.6],     2/√3),
        ("Rhombohedral (α=70°)",
            rhombohedral_basis(1.0, 70π/180),
            [1.0, 1.0, 1.0],     rhombohedral_defect(70π/180)),
    ]

    for (name, A, expected_norms, expected_defect) ∈ bravais
        # Sanity-check the table entry itself: if this fails, the expected
        # values don't match the primitive basis and the test is useless.
        @test sort(norm.(eachcol(A))) ≈ sort(expected_norms) rtol=1e-12
        @test orthogonalityDefect(A) ≈ expected_defect rtol=1e-12
        # Reducer must recover those invariants from any skewed basis for
        # the same lattice.
        for i ∈ 1:20
            M = A * RandUnimodMat3(4)
            R = minkReduce(M)
            @test isMinkReduced(R)
            @test sort(norm.(eachcol(R))) ≈ sort(expected_norms) rtol=1e-10
            @test orthogonalityDefect(R) ≈ expected_defect rtol=1e-10
        end
    end
end
