using MinkowskiReduction
using Test
using Random
using LinearAlgebra
using Documenter: DocMeta, doctest

# Run the jldoctests embedded in docstrings as part of the main test suite,
# so `Pkg.test()` validates them locally without a separate docs build.
@testset "Doctests" begin
    DocMeta.setdocmeta!(
        MinkowskiReduction,
        :DocTestSetup,
        :(using MinkowskiReduction; using LinearAlgebra);
        recursive=true,
    )
    doctest(MinkowskiReduction)
end

@testset "MinkowskiReduction.jl" begin
    U=[1, 2, 3];V=[-1, 2, 3];W=[3, 0, 4]
    # mink_reduce now returns (U, V, W, P, n); the original test stripped P out of the tuple.
    let (Ur, Vr, Wr, _, n) = mink_reduce(U, V, W)
        @test all((Ur, Vr, Wr, n) .≈ ([-2.0, 0.0, 0.0], [0.0, -2.0, 1.0], [-1.0, 2.0, 3.0], 2))
    end
    @test orthogonality_defect(mink_reduce(U,V,W)[1:3]...)≈1.0458250331675945
    @test orthogonality_defect(U,V,W)==4.375
    fcc = [[1,1,0],[1,0,1],[0,1,1]]
    @test orthogonality_defect(fcc...)≈1.4142135623730954
    @test orthogonality_defect(hcat(fcc...))≈1.4142135623730954
    @test orthogonality_defect(mink_reduce(U,V,W)[1:3]...)≈1.0458250331675945
    @test orthogonality_defect(hcat([U,V,W]...))==4.375
    m = devious_mat(26) # Largest size that doesn't overflow
    @test abs(det(hcat(mink_reduce(m[:,1],m[:,2],m[:,3])[1:3]...)))==1
    # shorten_w_in_uvw now returns (U, V, W, δP); original test pre-dates P.
    @test all(MinkowskiReduction.shorten_w_in_uvw(U,V,W)[1:3] .≈ ([-2.0, 0.0, 0.0], [1.0, 2.0, 3.0], [0.0, -2.0, 1.0]))
    # devious_mat stress tests. The reduced basis is a cubic cell (norms all
    # 1, |det|=1); the specific sign/permutation depends on rounding
    # details of Float64 arithmetic, so we assert only the lattice-intrinsic
    # invariants. `m * P == R` must hold exactly because P is exact integer.
    let m = devious_mat(26), (a, b, c, P, n) = mink_reduce(m[:,1], m[:,2], m[:,3])
        @test sort(norm.([a, b, c])) ≈ [1, 1, 1]
        @test is_mink_reduced(a, b, c)
        @test abs(det(BigInt.(P))) == 1
        @test m * P == hcat(a, b, c)
        @test 1 ≤ n ≤ 15
    end
    let m = devious_mat(20), (a, b, c, P, n) = mink_reduce(m[:,1], m[:,2], m[:,3])
        @test sort(norm.([a, b, c])) ≈ [1, 1, 1]
        @test is_mink_reduced(a, b, c)
        @test abs(det(BigInt.(P))) == 1
        @test m * P == hcat(a, b, c)
        @test 1 ≤ n ≤ 15
    end
    # These next few are not robust, but testing something that depends on the random number generator is difficult
    @test det(MinkowskiReduction.rand_lower_tri(35))≈1.0
    @test det(MinkowskiReduction.rand_unimod_mat2(4))≈1.0
    @test det(MinkowskiReduction.rand_lower_tri(45))≈1.0
    @test det(MinkowskiReduction.rand_unimod_mat2(6))≈1.0
    @test (MinkowskiReduction.fibonacci_mat(45)==[1836311903 2971215073; 1134903170 1836311903])
    @test !(U = [1, 2, 9]; V = [-3, 4, 3]; W = [3, 0, 4]; is_mink_reduced(U,V,W))
    @test !(U = [1, 2, 3]; V = [-3, 4, 3]; W = [3, 0, 4]; is_mink_reduced(U,V,W))
    @test !(U = [1, 0, 0]; V = [-1, 0, 1]; W = [3, 0, 4]; is_mink_reduced(U,V,W))
    @test !(U = [1, 2, 3]; V = [-1, 2, 3]; W = [3, 0, 4]; is_mink_reduced(U,V,W))
    @test !(U = [1, 0, 0]; V = [0,-4,-3]; W = [-3, 0, 4]; is_mink_reduced(U,V,W))
    @test !(U = [1, 0, 0]; V = [0, 4, 3]; W = [3, 0, 4]; is_mink_reduced(U,V,W))
    @test !(U = [0, 4, 3]; V = [-5, 0, 0]; W = [3, 0, 4]; is_mink_reduced(U,V,W))
    @test !(U = [1,0,0]; V = [0,1,0]; W = [0,1,-2]; is_mink_reduced(U,V,W))
    @test !(U = [1,0,-2]; V = [-1,1,-2]; W = [0,0,3]; is_mink_reduced(U,V,W))
    @test !(δ=.50;U = [1,0,0]; V = [-.4,1,0]; W = [-δ,-δ,5]; is_mink_reduced(U,V,W))
    @test !(U = [2,0,0]; V = [-.9,2,0]; W = [-1,-1,10]; is_mink_reduced(U,V,W))
    @test !(U = [2,0,0]; V = -[-.9,2,0]; W = [-1,-1,10]; is_mink_reduced(U,V,W))
    @test !(U = -[2,0,0]; V = -[-.9,2,0]; W = [-1,-1,10]; is_mink_reduced(U,V,W))
    @test !(U = -[2,0,0]; V = [-.9,2,0]; W = [-1,-1,10]; is_mink_reduced(U,V,W))
    @test (U = -[2,0,0]; V = [0,2,0]; W = [-1,-1,10]; is_mink_reduced(U,V,W))
    U = -[2,0,0]; V = [0,2,0]; W=U+V # Linearly dependent basis
    @test_throws ErrorException mink_reduce(U,V,W)
    @test_throws ErrorException mink_reduce(hcat(U,V,W))
    U = -[2.,0,0]; V = [0,2.,0]; W=U+V+[0,0,1e-170]
    @test_throws ErrorException mink_reduce(U,V,W)
    @test_throws ErrorException mink_reduce(hcat(U,V,W))
    @test (U = -[2.,0,0]; V = [0,2.,0]; W=U+V+[0,0,1e-150];
           let (a,b,c,_,n) = mink_reduce(U,V,W)
               (a,b,c,n) == ([0.0, 0.0, 1.0e-150], [-2.0, -0.0, -0.0], [0.0, 2.0, 0.0], 2)
           end)
    @test (U = -[2.,0,0]; V = [0,2.,0]; W=U+V+[0,0,1e-150];
           mink_reduce(hcat(U,V,W))[1] ==
           hcat([0.0, 0.0, 1.0e-150], [-2.0, -0.0, -0.0], [0.0, 2.0, 0.0]))
    #U = -[2^63-1,0,0]; V = [0,2^63-1,0]; W=U+V+[0,0,1e-150];
    #@test_throws InexactError mink_reduce(U,V,W) # silently overflows with later versions of julia
    for i ∈ 1:100
        @test is_permutation_matrix(vcat(shuffle!([[1 0 0],[0 1 0],[0 0 1]])...))
    end
    for i ∈ 1:100
        @test !is_permutation_matrix(vcat(shuffle!([[1 0 1],[0 1 0],[0 0 1]])...))
        @test !is_permutation_matrix(vcat(shuffle!([[1 0 0],[0 0 0],[0 0 1]])...))
    end
    A = [-1.7233692904465637e-9 0.0025000286163305314 0.0024999524070139123; 0.0024999730832105326 2.591951474986415e-8 0.0025000013059337757; 0.0024999629647447772 0.002499967442175358 -1.891891458670977e-8]
    for a ∈ logrange(1e-15,1e15,30)
        @test is_mink_reduced(mink_reduce(A*a)[1])
    end
    for i ∈ 1:10
        M_transform = rand_unimod_mat3(10)
        @test is_mink_reduced(mink_reduce(M_transform*A)[1])
    end
    A = [0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0]
    for ε ∈ logrange(1e-5,1e-1,10)
        for i ∈ 1:100
            noise = (2*rand(3,3).-1)*ε
            @test is_mink_reduced(mink_reduce(A+noise)[1])
            @test is_mink_reduced(mink_reduce(hcat(eachcol(A+noise)...))[1])
        end
    end
    # Test of noise levesls for BCC
    A = [1.0 1.0 -1.0; 1.0 1.0 1.0; -1.0 1.0 1.0]
    for ε ∈ logrange(1e-5,1.1e-1,20)
        for i ∈ 1:100
            noise = (2*rand(3,3).-1)*ε
            @test is_mink_reduced(mink_reduce(A+noise)[1])
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
                @test is_mink_reduced(mink_reduce(A+noise)[1])
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
                @test is_mink_reduced(mink_reduce(A+noise)[1])
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

    # --- 1. mink_reduce invariants on random unimodular transforms ------------
    # For each canonical Bravais lattice, apply random unimodular transforms
    # and check: the reduced output satisfies all 12 Minkowski conditions,
    # |det| is preserved, and the orthogonality defect does not increase.
    for A ∈ seed_bases
        for i ∈ 1:50
            # k=4 keeps entries of rand_unimod_mat3 small enough that `det` of
            # the resulting matrix is accurate to machine precision; larger
            # k would exercise the reducer equally well but introduce
            # spurious cancellation error in the reference det itself.
            M = A * rand_unimod_mat3(4)
            R, P = mink_reduce(M)
            @test is_mink_reduced(R)
            @test isapprox(abs(det(R)), abs(det(M)); rtol = 1e-10)
            @test orthogonality_defect(R) ≤ orthogonality_defect(M) * (1 + 1e-10)
            # Transform matrix invariants (added when P was threaded in):
            # R == M * P, P is integer and unimodular.
            @test eltype(P) <: Integer
            @test abs(det(BigInt.(P))) == 1   # exact: det(Float64(P)) loses precision when entries are large
            @test M * P ≈ R
        end
    end

    # --- 2. Idempotence ------------------------------------------------------
    # Reducing an already-reduced basis must be a no-op, and reducing it a
    # second time must yield the identity transform.
    for i ∈ 1:50
        B, _ = mink_reduce(Float64.(rand_unimod_mat3(8)))
        B2, P2 = mink_reduce(B)
        @test B2 ≈ B
        @test P2 == Matrix{Int}(I, 3, 3)
    end

    # --- 3. gauss_reduce (currently has no dedicated tests) -------------------
    # 2D Minkowski conditions on output, |det| preservation, idempotence,
    # error on linearly dependent input, and P-matrix invariants.
    for i ∈ 1:100
        U = rand(-100:100, 2)
        V = rand(-100:100, 2)
        (iszero(U) || iszero(V)) && continue
        d = abs(det(hcat(U, V)))
        d == 0 && continue                   # skip linearly dependent inputs
        short, long, P = gauss_reduce(U, V)   # returns (shorter, longer, 2×2 transform)
        @test norm(short) ≤ norm(long) + sqrt(eps())
        @test norm(long) ≤ norm(long + short) + sqrt(eps())
        @test norm(long) ≤ norm(long - short) + sqrt(eps())
        @test abs(det(hcat(short, long))) ≈ d
        @test eltype(P) <: Integer
        @test abs(det(BigInt.(P))) == 1
        @test hcat(U, V) * P ≈ hcat(short, long)
    end
    # Idempotence on the docstring example; a second call yields the identity
    # transform since the input is already reduced.
    let (a1, b1, _) = gauss_reduce([5, 8], [8, 13])
        a2, b2, P2 = gauss_reduce(a1, b1)
        @test a2 ≈ a1 && b2 ≈ b1
        @test P2 == Matrix{Int}(I, 2, 2)
    end
    # Parallel vectors are detected and raise an error
    @test_throws ErrorException gauss_reduce([1.0, 0, 0], [2.0, 0, 0])

    # --- 4. Sign and permutation invariance of the norm multiset -------------
    # Signed permutations of the input columns preserve the lattice, and so
    # must leave the sorted sequence of reduced norms unchanged (even though
    # the output column order / signs may differ due to non-uniqueness).
    norm_multiset(M) = sort(norm.(eachcol(M)))
    perms_3 = [[1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]]
    for A ∈ seed_bases
        ref = norm_multiset(mink_reduce(A)[1])
        for p ∈ perms_3, s1 ∈ (-1,1), s2 ∈ (-1,1), s3 ∈ (-1,1)
            Ap = A[:, p] .* reshape([s1, s2, s3], 1, 3)
            @test norm_multiset(mink_reduce(Ap)[1]) ≈ ref
        end
    end

    # --- 5. Scale invariance of the iteration count --------------------------
    # The algorithm makes only ratio-based decisions, so uniformly scaling
    # the basis must not change how many outer iterations it takes.
    U0, V0, W0 = [1.0, 2, 3], [-1.0, 2, 3], [3.0, 0, 4]
    _, _, _, _, n_ref = mink_reduce(U0, V0, W0)
    for α ∈ logrange(1e-8, 1e8, 20)
        _, _, _, _, n = mink_reduce(α*U0, α*V0, α*W0)
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
        @test orthogonality_defect(A) ≈ expected_defect rtol=1e-12
        # Reducer must recover those invariants from any skewed basis for
        # the same lattice.
        for i ∈ 1:20
            M = A * rand_unimod_mat3(4)
            R, P = mink_reduce(M)
            @test is_mink_reduced(R)
            @test sort(norm.(eachcol(R))) ≈ sort(expected_norms) rtol=1e-10
            @test orthogonality_defect(R) ≈ expected_defect rtol=1e-10
            @test M * P ≈ R
            @test abs(det(BigInt.(P))) == 1   # exact: det(Float64(P)) loses precision when entries are large
        end
    end

    # ------------------------------------------------------------------------
    # 7. Transform matrix P — dedicated sanity sweep.
    # The contract for downstream consumers (e.g. Spacey.jl) is:
    #   - `R == M * P` exactly (float equality up to Float64 roundoff),
    #   - `P` is an integer 3×3 matrix with `|det(P)| = 1`,
    #   - the three-vector form returns the same `P` as the matrix form.
    # ------------------------------------------------------------------------
    p_sweep_bases = [Matrix{Float64}(I, 3, 3),
                     [0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0],
                     [-0.5 0.5 0.5; 0.5 -0.5 0.5; 0.5 0.5 -0.5]]
    for i ∈ 1:100
        A = p_sweep_bases[rand(1:3)]
        M = A * rand_unimod_mat3(rand(1:8))
        R, P = mink_reduce(M)
        @test eltype(P) <: Integer
        @test size(P) == (3, 3)
        @test abs(det(BigInt.(P))) == 1
        @test M * P ≈ R
    end
    # Three-vector and matrix forms must agree on P.
    for i ∈ 1:20
        M = Float64.(rand_unimod_mat3(5))
        U, V, W, P_vec, _ = mink_reduce(M[:,1], M[:,2], M[:,3])
        R_mat, P_mat = mink_reduce(M)
        @test P_vec == P_mat
        @test hcat(U, V, W) == R_mat
    end

    # ------------------------------------------------------------------------
    # 8. Pure random Gaussian inputs.
    # All other random tests start from a structured Bravais lattice and apply
    # a `rand_unimod_mat3` (|det| = 1) transform — they exercise lattices with a
    # discrete, integer-valued underlying structure. Unstructured `randn`
    # inputs cover a different regime: arbitrary determinant, no integer
    # relationship between input and output, no a-priori bound on aspect
    # ratio. The reducer's contract still applies: output must satisfy the
    # 12 Minkowski conditions, |det| must be preserved exactly (the algorithm
    # only does integer linear combinations), and the integer transform `P`
    # must satisfy `M * P ≈ R` with `|det(P)| = 1`.
    #
    # We use a fixed RNG seed so the test is reproducible across runs but
    # samples a fresh distribution every release. Rare singular draws are
    # filtered by determinant magnitude; pathologically high aspect ratios
    # may hit mink_reduce's 15-iteration cap, in which case the call errors
    # — we count those and require the failure rate to be negligible.
    # ------------------------------------------------------------------------
    @testset "Pure random Gaussian inputs" begin
        rng = MersenneTwister(20260429)
        n_attempted = 100
        n_skipped_singular = 0
        n_iter_cap = 0
        n_passed = 0
        for _ ∈ 1:n_attempted
            M = randn(rng, 3, 3)
            if abs(det(M)) < 1e-6
                n_skipped_singular += 1
                continue
            end
            try
                R, P = mink_reduce(M)
                @test is_mink_reduced(R)
                @test isapprox(abs(det(R)), abs(det(M)); rtol = 1e-8)
                @test eltype(P) <: Integer
                @test abs(det(BigInt.(P))) == 1
                @test M * P ≈ R
                n_passed += 1
            catch e
                # The known failure mode is the 15-iteration cap; surface
                # anything else as a real bug.
                if e isa ErrorException && occursin("iter", e.msg)
                    n_iter_cap += 1
                else
                    rethrow(e)
                end
            end
        end
        # Almost all draws should pass. Singular and iter-cap failures are
        # both expected to be zero or very small. If the iter-cap rate
        # creeps up above ~5%, that's a regression worth investigating.
        @test n_passed ≥ Int(round(0.95 * n_attempted))
        @test n_iter_cap ≤ Int(round(0.05 * n_attempted))
    end
end
