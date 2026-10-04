using Test
using LinearAlgebra
using PfaPack

include("test_complex_rank2_c_order.jl")
include("test_real_rank2_c_order.jl")

@testset "PfaPack.jl" begin
    @testset "pfaffian_ltl!" begin
        @testset "Real matrices" begin
            # Test 2x2 skew-symmetric matrix
            A2 = [0.0 1.0; -1.0 0.0]
            A2_copy = copy(A2)
            pf2 = pfaffian_ltl!(A2_copy)
            @test pf2 ≈ 1.0

            # Test 4x4 skew-symmetric matrix
            A4 = zeros(4, 4)
            A4[1, 2] = 1.0
            A4[1, 3] = 2.0
            A4[1, 4] = 3.0
            A4[2, 3] = 4.0
            A4[2, 4] = 5.0
            A4[3, 4] = 6.0
            # Make it skew-symmetric
            for i in 1:4
                for j in (i+1):4
                    A4[j, i] = -A4[i, j]
                end
            end
            A4_copy = copy(A4)
            pf4 = pfaffian_ltl!(A4_copy)
            # Pfaffian of 4x4 matrix: pf(A) = A[1,2]*A[3,4] - A[1,3]*A[2,4] + A[1,4]*A[2,3]
            expected_pf4 = A4[1, 2]*A4[3, 4] - A4[1, 3]*A4[2, 4] + A4[1, 4]*A4[2, 3]
            @test pf4 ≈ expected_pf4

            # Test odd-sized matrix (should return 0)
            A3 = [0.0 1.0 2.0; -1.0 0.0 3.0; -2.0 -3.0 0.0]
            A3_copy = copy(A3)
            pf3 = pfaffian_ltl!(A3_copy)
            @test pf3 == 0.0
        end

        @testset "Complex matrices" begin
            # Test 2x2 complex skew-symmetric matrix
            A2c = [0.0+0.0im 1.0+1.0im; -1.0-1.0im 0.0+0.0im]
            A2c_copy = copy(A2c)
            pf2c = pfaffian_ltl!(A2c_copy)
            # For 2x2: pf(A) = A[1,2]
            @test pf2c ≈ A2c[1, 2]
        end

        @testset "overwrite_a parameter" begin
            # Use a larger matrix that will be modified
            A = zeros(4, 4)
            A[1, 2] = 1.0
            A[1, 3] = 2.0
            A[1, 4] = 3.0
            A[2, 3] = 4.0
            A[2, 4] = 5.0
            A[3, 4] = 6.0
            for i in 1:4
                for j in (i+1):4
                    A[j, i] = -A[i, j]
                end
            end
            A_orig = copy(A)
            pf1 = pfaffian_ltl!(A; overwrite_a=true)
            # Matrix should be modified for larger matrices
            @test A ≠ A_orig

            A2 = zeros(4, 4)
            A2[1, 2] = 1.0
            A2[1, 3] = 2.0
            A2[1, 4] = 3.0
            A2[2, 3] = 4.0
            A2[2, 4] = 5.0
            A2[3, 4] = 6.0
            for i in 1:4
                for j in (i+1):4
                    A2[j, i] = -A2[i, j]
                end
            end
            A2_orig = copy(A2)
            pf2 = pfaffian_ltl!(A2; overwrite_a=false)
            # Matrix should not be modified (internal copy is used)
            @test A2 == A2_orig
            @test pf1 ≈ pf2
        end
    end

    @testset "utu2pfa" begin
        @testset "Simple cases" begin
            # Test with identity permutation
            n = 4
            A = zeros(4, 4)
            # Upper triangular with unit diagonal (LTL decomposed form)
            A[1, 2] = 1.0
            A[1, 3] = 0.0
            A[1, 4] = 0.0
            A[2, 3] = 2.0
            A[2, 4] = 0.0
            A[3, 4] = 3.0
            iPiv = [1, 2, 3, 4]  # No permutation
            pf = utu2pfa(n, A, n, iPiv)
            # For utu2pfa: pf = A[1,2] * A[3,4] (for pairs)
            expected_pf = A[1, 2] * A[3, 4]
            @test pf ≈ expected_pf

            # Test with permutation
            iPiv2 = [2, 1, 4, 3]  # Swap 1<->2 and 3<->4
            pf2 = utu2pfa(n, A, n, iPiv2)
            # Each swap changes sign (2 swaps = sign change * 2 = no change)
            # Actually, each swap changes sign, so 2 swaps = no net change
            @test pf2 ≈ expected_pf
        end
    end

    @testset "utu2inv!" begin
        @testset "Simple 2x2 case" begin
            n = 2
            # Upper triangular LTL decomposed matrix
            A = [0.0 1.0; 0.0 0.0]  # Unit upper triangular
            A_orig = copy(A)
            iPiv = [1, 2]
            vT = zeros(n-1)
            M = zeros(n, n)

            utu2inv!(n, A, n, iPiv, vT, M, n)

            # For 2x2, the inverse should satisfy A * A_inv ≈ I (with skew-symmetry)
            # This is a basic test - full validation would require more complex setup
            @test size(A) == (n, n)
            @test size(M) == (n, n)
        end

        @testset "4x4 case" begin
            n = 4
            A = zeros(n, n)
            # Upper triangular part (unit diagonal)
            A[1, 2] = 1.0
            A[1, 3] = 0.5
            A[1, 4] = 0.3
            A[2, 3] = 2.0
            A[2, 4] = 0.7
            A[3, 4] = 3.0
            A_orig = copy(A)
            iPiv = [1, 2, 3, 4]
            vT = zeros(n-1)
            M = zeros(n, n)

            utu2inv!(n, A, n, iPiv, vT, M, n)

            # Basic size checks
            @test size(A) == (n, n)
            @test size(M) == (n, n)
            @test length(vT) == n-1
        end
    end

    @testset "utu2inv! vs cimpl_utu2inv! comparison" begin
        @testset "Float64 matrices" begin
            for n in [4, 6, 8, 10]
                # Create random upper triangular LTL decomposed matrix
                A_julia = zeros(Float64, n, n)
                for i in 1:n
                    for j in (i+1):n
                        A_julia[i, j] = randn()
                    end
                end
                A_c = copy(A_julia)

                iPiv = collect(1:n)
                vT_julia = zeros(Float64, n-1)
                vT_c = zeros(Float64, n-1)
                M_julia = zeros(Float64, n, n)
                M_c = zeros(Float64, n, n)

                # Call Julia implementation
                utu2inv!(n, A_julia, n, iPiv, vT_julia, M_julia, n)

                # Call C++ implementation
                cimpl_utu2inv!(n, A_c, n, iPiv, vT_c, M_c, n)

                # Compare results
                @test A_julia ≈ A_c rtol=1e-12 atol=1e-14
                @test M_julia ≈ M_c rtol=1e-12 atol=1e-14
                @test vT_julia ≈ vT_c rtol=1e-12 atol=1e-14
            end
        end

        @testset "Complex{Float64} matrices" begin
            for n in [4, 6, 8, 10]
                # Create random upper triangular LTL decomposed matrix
                A_julia = zeros(Complex{Float64}, n, n)
                for i in 1:n
                    for j in (i+1):n
                        A_julia[i, j] = randn() + im * randn()
                    end
                end
                A_c = copy(A_julia)

                iPiv = collect(1:n)
                vT_julia = zeros(Complex{Float64}, n-1)
                vT_c = zeros(Complex{Float64}, n-1)
                M_julia = zeros(Complex{Float64}, n, n)
                M_c = zeros(Complex{Float64}, n, n)

                # Call Julia implementation
                utu2inv!(n, A_julia, n, iPiv, vT_julia, M_julia, n)

                # Call C++ implementation
                cimpl_utu2inv!(n, A_c, n, iPiv, vT_c, M_c, n)

                # Compare results
                @test A_julia ≈ A_c rtol=1e-12 atol=1e-14
                @test M_julia ≈ M_c rtol=1e-12 atol=1e-14
                @test vT_julia ≈ vT_c rtol=1e-12 atol=1e-14
            end
        end

        @testset "With non-trivial pivot" begin
            n = 6
            A_julia = zeros(Complex{Float64}, n, n)
            for i in 1:n
                for j in (i+1):n
                    A_julia[i, j] = randn() + im * randn()
                end
            end
            A_c = copy(A_julia)

            # Non-trivial pivot (swaps within pairs)
            iPiv = [2, 1, 4, 3, 6, 5]
            vT_julia = zeros(Complex{Float64}, n-1)
            vT_c = zeros(Complex{Float64}, n-1)
            M_julia = zeros(Complex{Float64}, n, n)
            M_c = zeros(Complex{Float64}, n, n)

            # Call Julia implementation
            utu2inv!(n, A_julia, n, iPiv, vT_julia, M_julia, n)

            # Call C++ implementation
            cimpl_utu2inv!(n, A_c, n, iPiv, vT_c, M_c, n)

            # Compare results
            @test A_julia ≈ A_c rtol=1e-12 atol=1e-14
            @test M_julia ≈ M_c rtol=1e-12 atol=1e-14
            @test vT_julia ≈ vT_c rtol=1e-12 atol=1e-14
        end

        @testset "With actual LTL decomposition result" begin
            # Use actual LTL decomposition output as input
            n = 6
            A_orig = zeros(Complex{Float64}, n, n)
            for i in 1:n
                for j in (i+1):n
                    A_orig[i, j] = randn() + im * randn()
                    A_orig[j, i] = -A_orig[i, j]  # Skew-symmetric
                end
            end

            A_ltl = copy(A_orig)
            iPiv = zeros(Int, n)
            info = julia_zsktf2!(A_ltl, iPiv)

            if info == 0  # Decomposition succeeded
                A_julia = copy(A_ltl)
                A_c = copy(A_ltl)
                vT_julia = zeros(Complex{Float64}, n-1)
                vT_c = zeros(Complex{Float64}, n-1)
                M_julia = zeros(Complex{Float64}, n, n)
                M_c = zeros(Complex{Float64}, n, n)

                # Call Julia implementation
                utu2inv!(n, A_julia, n, iPiv, vT_julia, M_julia, n)

                # Call C++ implementation
                cimpl_utu2inv!(n, A_c, n, iPiv, vT_c, M_c, n)

                # Compare results
                @test A_julia ≈ A_c rtol=1e-12 atol=1e-14
                @test M_julia ≈ M_c rtol=1e-12 atol=1e-14
                @test vT_julia ≈ vT_c rtol=1e-12 atol=1e-14
            end
        end
    end

    @testset "Integration tests" begin
        @testset "pfaffian_ltl! consistency" begin
            # Test that pfaffian_ltl! gives consistent results
            A = zeros(4, 4)
            A[1, 2] = 1.0
            A[1, 3] = 2.0
            A[1, 4] = 3.0
            A[2, 3] = 4.0
            A[2, 4] = 5.0
            A[3, 4] = 6.0
            # Make skew-symmetric
            for i in 1:4
                for j in (i+1):4
                    A[j, i] = -A[i, j]
                end
            end

            A1 = copy(A)
            A2 = copy(A)
            pf1 = pfaffian_ltl!(A1, overwrite_a=true)
            pf2 = pfaffian_ltl!(A2, overwrite_a=true)

            # Should give same result for same input
            @test pf1 ≈ pf2
        end
    end

    @testset "Edge cases" begin
        @testset "Zero matrix" begin
            A = zeros(4, 4)
            pf = pfaffian_ltl!(A)
            @test pf == 0.0
        end

        @testset "Non-square matrix" begin
            A = zeros(3, 4)
            @test_throws ArgumentError pfaffian_ltl!(A)
        end
    end

    @testset "julia_zsktf2! vs fimpl_zsktf2! comparison" begin
        @testset "Real matrices (converted to Complex)" begin
            # Test 4x4 skew-symmetric matrix
            # Note: fimpl_zsktf2! is for complex matrices, so we use ComplexF64
            n = 4
            A_orig = zeros(Float64, n, n)
            A_orig[1, 2] = 1.0
            A_orig[1, 3] = 2.0
            A_orig[1, 4] = 3.0
            A_orig[2, 3] = 4.0
            A_orig[2, 4] = 5.0
            A_orig[3, 4] = 6.0
            # Make it skew-symmetric
            for i in 1:n
                for j in (i+1):n
                    A_orig[j, i] = -A_orig[i, j]
                end
            end

            A_julia = ComplexF64.(A_orig)
            A_fortran = copy(A_julia)
            iPiv_julia = zeros(Int, n)
            iPiv_fortran = zeros(Cint, n)

            # Call Julia implementation
            info_julia = julia_zsktf2!(A_julia, iPiv_julia)

            # Call Fortran implementation (mode='N' for normal mode, matching julia_zsktf2!)
            info_fortran = fimpl_zsktf2!(A_fortran, iPiv_fortran; uplo='U', mode='N')

            # Both should return same info
            @test info_julia == info_fortran

            # Compare upper triangular part only (UPLO='U' means only upper part is computed)
            # Lower triangular part is not referenced by Fortran implementation
            for i in 1:n
                for j in i:n  # Only upper triangular (including diagonal)
                    @test A_julia[i, j] ≈ A_fortran[i, j] rtol=1e-12 atol=1e-14
                end
            end

            # Compare pivot arrays (convert Cint to Int for comparison)
            @test iPiv_julia == Int.(iPiv_fortran)
        end

        @testset "Complex matrices" begin
            # Test 4x4 complex skew-symmetric matrix
            n = 4
            A_julia = zeros(ComplexF64, n, n)
            A_julia[1, 2] = 1.0 + 0.5im
            A_julia[1, 3] = 2.0 + 1.0im
            A_julia[1, 4] = 3.0 + 1.5im
            A_julia[2, 3] = 4.0 + 2.0im
            A_julia[2, 4] = 5.0 + 2.5im
            A_julia[3, 4] = 6.0 + 3.0im
            # Make it skew-symmetric
            for i in 1:n
                for j in (i+1):n
                    A_julia[j, i] = -A_julia[i, j]
                end
            end

            A_fortran = copy(A_julia)
            iPiv_julia = zeros(Int, n)
            iPiv_fortran = zeros(Cint, n)

            # Call Julia implementation
            info_julia = julia_zsktf2!(A_julia, iPiv_julia)

            # Call Fortran implementation (mode='N' for normal mode)
            info_fortran = fimpl_zsktf2!(A_fortran, iPiv_fortran; uplo='U', mode='N')

            # Both should return same info
            @test info_julia == info_fortran

            # Compare upper triangular part only
            for i in 1:n
                for j in i:n
                    @test A_julia[i, j] ≈ A_fortran[i, j] rtol=1e-12 atol=1e-14
                end
            end

            # Compare pivot arrays
            @test iPiv_julia == Int.(iPiv_fortran)
        end

        @testset "Larger matrices" begin
            # Test 6x6 skew-symmetric matrix
            n = 6
            A_julia = zeros(ComplexF64, n, n)
            # Fill upper triangular part with random-like values
            A_julia[1, 2] = 1.0 + 0.1im
            A_julia[1, 3] = 2.0 + 0.2im
            A_julia[1, 4] = 3.0 + 0.3im
            A_julia[1, 5] = 4.0 + 0.4im
            A_julia[1, 6] = 5.0 + 0.5im
            A_julia[2, 3] = 6.0 + 0.6im
            A_julia[2, 4] = 7.0 + 0.7im
            A_julia[2, 5] = 8.0 + 0.8im
            A_julia[2, 6] = 9.0 + 0.9im
            A_julia[3, 4] = 10.0 + 1.0im
            A_julia[3, 5] = 11.0 + 1.1im
            A_julia[3, 6] = 12.0 + 1.2im
            A_julia[4, 5] = 13.0 + 1.3im
            A_julia[4, 6] = 14.0 + 1.4im
            A_julia[5, 6] = 15.0 + 1.5im
            # Make it skew-symmetric
            for i in 1:n
                for j in (i+1):n
                    A_julia[j, i] = -A_julia[i, j]
                end
            end

            A_fortran = copy(A_julia)
            iPiv_julia = zeros(Int, n)
            iPiv_fortran = zeros(Cint, n)

            # Call Julia implementation
            info_julia = julia_zsktf2!(A_julia, iPiv_julia)

            # Call Fortran implementation (mode='N' for normal mode)
            info_fortran = fimpl_zsktf2!(A_fortran, iPiv_fortran; uplo='U', mode='N')

            # Both should return same info
            @test info_julia == info_fortran

            # Compare upper triangular part only
            for i in 1:n
                for j in i:n
                    @test A_julia[i, j] ≈ A_fortran[i, j] rtol=1e-12 atol=1e-14
                end
            end

            # Compare pivot arrays
            @test iPiv_julia == Int.(iPiv_fortran)
        end

        @testset "Edge case: zero column" begin
            # Test matrix with a zero column (should set info > 0)
            n = 4
            A_julia = zeros(ComplexF64, n, n)
            A_julia[1, 2] = 1.0
            A_julia[1, 3] = 2.0
            # Column 4 is zero (except diagonal which is always 0 for skew-symmetric)
            A_julia[2, 3] = 3.0
            # Make it skew-symmetric
            for i in 1:n
                for j in (i+1):n
                    A_julia[j, i] = -A_julia[i, j]
                end
            end

            A_fortran = copy(A_julia)
            iPiv_julia = zeros(Int, n)
            iPiv_fortran = zeros(Cint, n)

            # Call Julia implementation
            info_julia = julia_zsktf2!(A_julia, iPiv_julia)

            # Call Fortran implementation (mode='N' for normal mode)
            info_fortran = fimpl_zsktf2!(A_fortran, iPiv_fortran; uplo='U', mode='N')

            # Both should return same info (non-zero if zero column detected)
            @test info_julia == info_fortran

            # Compare pivot arrays
            @test iPiv_julia == Int.(iPiv_fortran)
        end
    end

    @testset "julia_zsktf2_turbo! (LoopVectorization optimized)" begin
        @testset "Correctness vs julia_zsktf2!" begin
            for n in [4, 8, 16, 32, 64]
                A_orig = zeros(ComplexF64, n, n)
                for i in 1:n
                    for j in (i+1):n
                        A_orig[i, j] = randn() + im * randn()
                        A_orig[j, i] = -A_orig[i, j]
                    end
                end

                A_julia = copy(A_orig)
                A_turbo = copy(A_orig)
                iPiv_julia = zeros(Int, n)
                iPiv_turbo = zeros(Int, n)

                info_julia = julia_zsktf2!(A_julia, iPiv_julia)
                info_turbo = julia_zsktf2_turbo!(A_turbo, iPiv_turbo)

                @test info_julia == info_turbo
                @test iPiv_julia == iPiv_turbo

                # Compare upper triangular part
                for i in 1:n
                    for j in i:n
                        @test A_julia[i, j] ≈ A_turbo[i, j] rtol=1e-12 atol=1e-14
                    end
                end
            end
        end

        @testset "Correctness vs Fortran" begin
            for n in [4, 8, 16, 32, 64]
                A_orig = zeros(ComplexF64, n, n)
                for i in 1:n
                    for j in (i+1):n
                        A_orig[i, j] = randn() + im * randn()
                        A_orig[j, i] = -A_orig[i, j]
                    end
                end

                A_turbo = copy(A_orig)
                A_fortran = copy(A_orig)
                iPiv_turbo = zeros(Int, n)
                iPiv_fortran = zeros(Cint, n)

                info_turbo = julia_zsktf2_turbo!(A_turbo, iPiv_turbo)
                info_fortran = fimpl_zsktf2!(A_fortran, iPiv_fortran; uplo='U', mode='N')

                @test info_turbo == info_fortran
                @test iPiv_turbo == Int.(iPiv_fortran)

                # Compare upper triangular part
                for i in 1:n
                    for j in i:n
                        @test A_turbo[i, j] ≈ A_fortran[i, j] rtol=1e-10 atol=1e-12
                    end
                end
            end
        end

        @testset "Edge case: zero column" begin
            n = 4
            A_turbo = zeros(ComplexF64, n, n)
            A_turbo[1, 2] = 1.0
            A_turbo[1, 3] = 2.0
            A_turbo[2, 3] = 3.0
            for i in 1:n
                for j in (i+1):n
                    A_turbo[j, i] = -A_turbo[i, j]
                end
            end

            A_julia = copy(A_turbo)
            iPiv_julia = zeros(Int, n)
            iPiv_turbo = zeros(Int, n)

            info_julia = julia_zsktf2!(A_julia, iPiv_julia)
            info_turbo = julia_zsktf2_turbo!(A_turbo, iPiv_turbo)

            @test info_julia == info_turbo
            @test iPiv_julia == iPiv_turbo
        end

        @testset "Non-square matrix throws error" begin
            A = zeros(ComplexF64, 3, 4)
            iPiv = zeros(Int, 3)
            @test_throws ArgumentError julia_zsktf2_turbo!(A, iPiv)
        end
    end
end

