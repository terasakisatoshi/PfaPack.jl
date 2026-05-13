"""
Benchmark script for comparing julia_zsktf2! (Julia) vs fimpl_zsktf2! (Fortran)

Usage:
    julia --project=. benchmark/benchmark_ltl.jl [--complex] [sizes...]

Examples:
    julia --project=. benchmark/benchmark_ltl.jl              # Real matrices, default sizes
    julia --project=. benchmark/benchmark_ltl.jl 4 8 16 32    # Real matrices, custom sizes
    julia --project=. benchmark/benchmark_ltl.jl --complex    # Complex matrices, compare with Fortran
    julia --project=. benchmark/benchmark_ltl.jl --complex 4 8 16 32

Requirements:
    - BenchmarkTools.jl (add with `Pkg.add("BenchmarkTools")`)
"""

using Pkg
Pkg.activate(dirname(@__DIR__))

using BenchmarkTools
using LinearAlgebra
using Printf
using PfaPack
using PfaPack: fimpl_zsktf2!, fimpl_dsktf2!
using PfaPack.LTLDecomposition: julia_zsktf2!, julia_dsktf2!, julia_zsktf2_turbo!

# Configure BenchmarkTools
BenchmarkTools.DEFAULT_PARAMETERS.samples = 100
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5

"""
Create a random skew-symmetric matrix of size n×n
"""
function create_skew_symmetric(::Type{T}, n::Int) where T
    A = randn(T, n, n)
    # Make skew-symmetric: A = (A - A') / 2
    A = (A - transpose(A)) / 2
    return A
end

"""
Run benchmark for real matrices (Julia vs Fortran comparison)
"""
function benchmark_real(n::Int; verbose::Bool=true)
    if verbose
        println("\n" * "="^60)
        println("Matrix size: $n × $n (Float64)")
        println("="^60)
    end

    A_orig = create_skew_symmetric(Float64, n)
    A_julia = similar(A_orig)
    A_fortran = similar(A_orig)
    iPiv_julia = zeros(Int, n)
    iPiv_fortran = zeros(Cint, n)

    # Verify results match
    A_julia .= A_orig
    A_fortran .= A_orig
    fill!(iPiv_julia, 0)
    fill!(iPiv_fortran, 0)

    julia_dsktf2!(A_julia, iPiv_julia)
    fimpl_dsktf2!(A_fortran, iPiv_fortran; uplo='U', mode='N')

    upper_match = all(isapprox(A_julia[i,j], A_fortran[i,j], rtol=1e-10, atol=1e-12)
                      for i in 1:n for j in i:n)
    pivot_match = iPiv_julia == Int.(iPiv_fortran)

    if verbose
        println("Results match (upper triangular): $upper_match")
        println("Pivots match: $pivot_match")
        println()
    end

    # Benchmark Julia (optimized julia_dsktf2!)
    if verbose
        println("Benchmarking Julia (julia_dsktf2!)...")
    end
    bench_julia = @benchmark begin
        $A_julia .= $A_orig
        fill!($iPiv_julia, 0)
        julia_dsktf2!($A_julia, $iPiv_julia)
    end

    # Benchmark Fortran mode='N'
    if verbose
        println("Benchmarking Fortran (fimpl_dsktf2!, mode='N')...")
    end
    bench_fortran = @benchmark begin
        $A_fortran .= $A_orig
        fill!($iPiv_fortran, 0)
        fimpl_dsktf2!($A_fortran, $iPiv_fortran; uplo='U', mode='N')
    end

    return (n=n, julia=bench_julia, fortran=bench_fortran,
            upper_match=upper_match, pivot_match=pivot_match)
end

"""
Run benchmark for complex matrices (Julia, Turbo, and Fortran)
"""
function benchmark_complex(n::Int; verbose::Bool=true)
    if verbose
        println("\n" * "="^60)
        println("Matrix size: $n × $n (ComplexF64)")
        println("="^60)
    end

    A_orig = create_skew_symmetric(ComplexF64, n)
    A_julia = similar(A_orig)
    A_turbo = similar(A_orig)
    A_fortran = similar(A_orig)
    iPiv_julia = zeros(Int, n)
    iPiv_turbo = zeros(Int, n)
    iPiv_fortran = zeros(Cint, n)

    # Verify results match
    A_julia .= A_orig
    A_turbo .= A_orig
    A_fortran .= A_orig
    fill!(iPiv_julia, 0)
    fill!(iPiv_turbo, 0)
    fill!(iPiv_fortran, 0)

    julia_zsktf2!(A_julia, iPiv_julia)
    julia_zsktf2_turbo!(A_turbo, iPiv_turbo)
    fimpl_zsktf2!(A_fortran, iPiv_fortran; uplo='U', mode='N')

    upper_match_julia = all(isapprox(A_julia[i,j], A_fortran[i,j], rtol=1e-10, atol=1e-12)
                           for i in 1:n for j in i:n)
    upper_match_turbo = all(isapprox(A_turbo[i,j], A_fortran[i,j], rtol=1e-10, atol=1e-12)
                           for i in 1:n for j in i:n)
    pivot_match_julia = iPiv_julia == Int.(iPiv_fortran)
    pivot_match_turbo = iPiv_turbo == Int.(iPiv_fortran)

    if verbose
        println("Julia vs Fortran (upper triangular): $upper_match_julia")
        println("Turbo vs Fortran (upper triangular): $upper_match_turbo")
        println("Pivots match (Julia): $pivot_match_julia")
        println("Pivots match (Turbo): $pivot_match_turbo")
        println()
    end

    # Benchmark Julia
    if verbose
        println("Benchmarking Julia (julia_zsktf2!)...")
    end
    bench_julia = @benchmark begin
        $A_julia .= $A_orig
        fill!($iPiv_julia, 0)
        julia_zsktf2!($A_julia, $iPiv_julia)
    end

    # Benchmark Turbo
    if verbose
        println("Benchmarking Turbo (julia_zsktf2_turbo!)...")
    end
    bench_turbo = @benchmark begin
        $A_turbo .= $A_orig
        fill!($iPiv_turbo, 0)
        julia_zsktf2_turbo!($A_turbo, $iPiv_turbo)
    end

    # Benchmark Fortran mode='N'
    if verbose
        println("Benchmarking Fortran (fimpl_zsktf2!, mode='N')...")
    end
    bench_fortran = @benchmark begin
        $A_fortran .= $A_orig
        fill!($iPiv_fortran, 0)
        fimpl_zsktf2!($A_fortran, $iPiv_fortran; uplo='U', mode='N')
    end

    return (n=n, julia=bench_julia, turbo=bench_turbo, fortran=bench_fortran,
            upper_match_julia=upper_match_julia, upper_match_turbo=upper_match_turbo,
            pivot_match_julia=pivot_match_julia, pivot_match_turbo=pivot_match_turbo)
end

"""
Run benchmarks for real matrices (with Fortran comparison)
"""
function run_real_benchmarks(sizes::Vector{Int})
    println("="^60)
    println("LTL Decomposition Benchmark (Real Matrices)")
    println("Comparing: julia_zsktf2! (Julia) vs fimpl_dsktf2! (Fortran)")
    println("="^60)
    println("\nMatrix sizes: $sizes")
    println("Element type: Float64")

    results = [benchmark_real(n) for n in sizes]

    # Summary
    println("\n\n" * "="^60)
    println("SUMMARY TABLE")
    println("="^60)
    @printf("\n%-6s | %-12s | %-12s | %-8s | %-8s\n",
            "Size", "Julia (ms)", "Fortran (ms)", "F/J", "Match")
    println("-"^55)

    for r in results
        t_j = median(r.julia).time / 1e6
        t_f = median(r.fortran).time / 1e6
        match = r.upper_match && r.pivot_match
        @printf("%-6d | %-12.6f | %-12.6f | %-8.2f | %-8s\n",
                r.n, t_j, t_f, t_f/t_j, match ? "✓" : "✗")
    end

    println("\nLegend:")
    println("  Julia   = julia_dsktf2! (optimized Julia, upper triangular only)")
    println("  Fortran = fimpl_dsktf2! mode='N' (full tridiagonalization)")
    println("  F/J     = Fortran time / Julia time (>1 = Julia faster)")
    println("  Match   = Results match within tolerance")

    return results
end

"""
Run benchmarks for complex matrices (Julia, Turbo, and Fortran comparison)
"""
function run_complex_benchmarks(sizes::Vector{Int})
    println("="^60)
    println("LTL Decomposition Benchmark (Complex Matrices)")
    println("Comparing: Julia vs Turbo vs Fortran")
    println("="^60)
    println("\nMatrix sizes: $sizes")
    println("Element type: ComplexF64")

    results = [benchmark_complex(n) for n in sizes]

    # Summary
    println("\n\n" * "="^60)
    println("SUMMARY TABLE")
    println("="^60)
    @printf("\n%-6s | %-12s | %-12s | %-12s | %-8s | %-8s | %-6s | %-6s\n",
            "Size", "Julia (ms)", "Turbo (ms)", "Fortran (ms)", "F/J", "F/T", "J_ok", "T_ok")
    println("-"^90)

    for r in results
        t_j = median(r.julia).time / 1e6
        t_t = median(r.turbo).time / 1e6
        t_f = median(r.fortran).time / 1e6
        match_j = r.upper_match_julia && r.pivot_match_julia
        match_t = r.upper_match_turbo && r.pivot_match_turbo
        @printf("%-6d | %-12.6f | %-12.6f | %-12.6f | %-8.2f | %-8.2f | %-6s | %-6s\n",
                r.n, t_j, t_t, t_f, t_f/t_j, t_f/t_t, match_j ? "✓" : "✗", match_t ? "✓" : "✗")
    end

    println("\nLegend:")
    println("  Julia   = julia_zsktf2! (pure Julia)")
    println("  Turbo   = julia_zsktf2_turbo! (LoopVectorization.jl)")
    println("  Fortran = fimpl_zsktf2! mode='N' (full tridiagonalization)")
    println("  F/J     = Fortran time / Julia time (>1 = Julia faster)")
    println("  F/T     = Fortran time / Turbo time (>1 = Turbo faster)")
    println("  J_ok    = Julia results match Fortran")
    println("  T_ok    = Turbo results match Fortran")

    return results
end

# Main
if abspath(PROGRAM_FILE) == @__FILE__
    use_complex = "--complex" in ARGS
    sizes_args = filter(x -> x != "--complex", ARGS)
    sizes = isempty(sizes_args) ? [4, 8, 16, 32, 64, 128] : parse.(Int, sizes_args)

    if use_complex
        run_complex_benchmarks(sizes)
    else
        run_real_benchmarks(sizes)
    end
end
