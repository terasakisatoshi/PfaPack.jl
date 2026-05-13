"""
Benchmark script for comparing utu2inv! (Julia) vs cimpl_utu2inv! (C++)

Usage:
    julia --project=. benchmark/benchmark_utu2inv.jl [--real] [sizes...]

Examples:
    julia --project=. benchmark/benchmark_utu2inv.jl              # Complex matrices, default sizes
    julia --project=. benchmark/benchmark_utu2inv.jl 4 8 16 32    # Complex matrices, custom sizes
    julia --project=. benchmark/benchmark_utu2inv.jl --real       # Real matrices
    julia --project=. benchmark/benchmark_utu2inv.jl --real 4 8 16 32

Requirements:
    - BenchmarkTools.jl (add with `Pkg.add("BenchmarkTools")`)
"""

using Pkg
Pkg.activate(dirname(@__DIR__))

using BenchmarkTools
using LinearAlgebra
using Printf
using PfaPack
using PfaPack: utu2inv!, cimpl_utu2inv!
using PfaPack.LTLDecomposition: julia_zsktf2!

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
Prepare LTL decomposed matrix for utu2inv! benchmark
"""
function prepare_ltl_input(::Type{T}, n::Int) where T
    A_orig = create_skew_symmetric(T, n)
    A_ltl = copy(A_orig)
    iPiv = zeros(Int, n)
    info = julia_zsktf2!(A_ltl, iPiv)
    if info != 0
        error("LTL decomposition failed with info=$info")
    end
    return A_ltl, iPiv
end

"""
Run benchmark for a specific element type
"""
function benchmark_utu2inv(::Type{T}, n::Int; verbose::Bool=true) where T
    type_name = T <: Complex ? "ComplexF64" : "Float64"

    if verbose
        println("\n" * "="^60)
        println("Matrix size: $n × $n ($type_name)")
        println("="^60)
    end

    # Prepare input
    A_ltl_orig, iPiv = prepare_ltl_input(T, n)

    # Allocate work arrays
    A_julia = similar(A_ltl_orig)
    A_cpp = similar(A_ltl_orig)
    vT_julia = zeros(T, n-1)
    vT_cpp = zeros(T, n-1)
    M_julia = zeros(T, n, n)
    M_cpp = zeros(T, n, n)

    # Verify results match
    A_julia .= A_ltl_orig
    A_cpp .= A_ltl_orig
    fill!(vT_julia, zero(T))
    fill!(vT_cpp, zero(T))
    fill!(M_julia, zero(T))
    fill!(M_cpp, zero(T))

    utu2inv!(n, A_julia, n, iPiv, vT_julia, M_julia, n)
    cimpl_utu2inv!(n, A_cpp, n, iPiv, vT_cpp, M_cpp, n)

    A_match = isapprox(A_julia, A_cpp, rtol=1e-12, atol=1e-14)
    M_match = isapprox(M_julia, M_cpp, rtol=1e-12, atol=1e-14)
    vT_match = isapprox(vT_julia, vT_cpp, rtol=1e-12, atol=1e-14)

    if verbose
        println("Results match:")
        println("  A matrices: $A_match")
        println("  M matrices: $M_match")
        println("  vT vectors: $vT_match")
        println()
    end

    # Benchmark Julia implementation
    if verbose
        println("Benchmarking Julia (utu2inv!)...")
    end
    bench_julia = @benchmark begin
        $A_julia .= $A_ltl_orig
        fill!($vT_julia, zero($T))
        fill!($M_julia, zero($T))
        utu2inv!($n, $A_julia, $n, $iPiv, $vT_julia, $M_julia, $n)
    end

    # Benchmark C++ implementation
    if verbose
        println("Benchmarking C++ (cimpl_utu2inv!)...")
    end
    bench_cpp = @benchmark begin
        $A_cpp .= $A_ltl_orig
        fill!($vT_cpp, zero($T))
        fill!($M_cpp, zero($T))
        cimpl_utu2inv!($n, $A_cpp, $n, $iPiv, $vT_cpp, $M_cpp, $n)
    end

    return (n=n, julia=bench_julia, cpp=bench_cpp,
            A_match=A_match, M_match=M_match, vT_match=vT_match)
end

"""
Run benchmarks for complex matrices
"""
function run_complex_benchmarks(sizes::Vector{Int})
    println("="^70)
    println("utu2inv! Benchmark (Complex Matrices)")
    println("Comparing: utu2inv! (Julia) vs cimpl_utu2inv! (C++)")
    println("="^70)
    println("\nMatrix sizes: $sizes")
    println("Element type: ComplexF64")

    results = [benchmark_utu2inv(ComplexF64, n) for n in sizes]

    # Summary
    println("\n\n" * "="^70)
    println("SUMMARY TABLE")
    println("="^70)
    @printf("\n%-8s | %-14s | %-14s | %-10s | %-8s\n",
            "Size", "Julia (ms)", "C++ (ms)", "C++/J", "Match")
    println("-"^70)

    for r in results
        t_julia = median(r.julia).time / 1e6
        t_cpp = median(r.cpp).time / 1e6
        ratio = t_cpp / t_julia
        match = r.A_match && r.M_match && r.vT_match
        @printf("%-8d | %-14.6f | %-14.6f | %-10.2f | %-8s\n",
                r.n, t_julia, t_cpp, ratio, match ? "✓" : "✗")
    end

    println("\nLegend:")
    println("  Julia   = utu2inv! (pure Julia implementation)")
    println("  C++     = cimpl_utu2inv! (C++ via ccall)")
    println("  C++/J   = C++ time / Julia time (>1 = Julia faster)")
    println("  Match   = Results match within tolerance")

    # Memory allocation comparison
    println("\n" * "="^70)
    println("MEMORY ALLOCATION")
    println("="^70)
    @printf("\n%-8s | %-16s | %-16s | %-16s | %-16s\n",
            "Size", "Julia allocs", "Julia mem", "C++ allocs", "C++ mem")
    println("-"^80)

    for r in results
        @printf("%-8d | %-16d | %-16s | %-16d | %-16s\n",
                r.n,
                r.julia.allocs, Base.format_bytes(r.julia.memory),
                r.cpp.allocs, Base.format_bytes(r.cpp.memory))
    end

    return results
end

"""
Run benchmarks for real matrices
"""
function run_real_benchmarks(sizes::Vector{Int})
    println("="^70)
    println("utu2inv! Benchmark (Real Matrices)")
    println("Comparing: utu2inv! (Julia) vs cimpl_utu2inv! (C++)")
    println("="^70)
    println("\nMatrix sizes: $sizes")
    println("Element type: Float64")

    results = [benchmark_utu2inv(Float64, n) for n in sizes]

    # Summary
    println("\n\n" * "="^70)
    println("SUMMARY TABLE")
    println("="^70)
    @printf("\n%-8s | %-14s | %-14s | %-10s | %-8s\n",
            "Size", "Julia (ms)", "C++ (ms)", "C++/J", "Match")
    println("-"^70)

    for r in results
        t_julia = median(r.julia).time / 1e6
        t_cpp = median(r.cpp).time / 1e6
        ratio = t_cpp / t_julia
        match = r.A_match && r.M_match && r.vT_match
        @printf("%-8d | %-14.6f | %-14.6f | %-10.2f | %-8s\n",
                r.n, t_julia, t_cpp, ratio, match ? "✓" : "✗")
    end

    println("\nLegend:")
    println("  Julia   = utu2inv! (pure Julia implementation)")
    println("  C++     = cimpl_utu2inv! (C++ via ccall)")
    println("  C++/J   = C++ time / Julia time (>1 = Julia faster)")
    println("  Match   = Results match within tolerance")

    # Memory allocation comparison
    println("\n" * "="^70)
    println("MEMORY ALLOCATION")
    println("="^70)
    @printf("\n%-8s | %-16s | %-16s | %-16s | %-16s\n",
            "Size", "Julia allocs", "Julia mem", "C++ allocs", "C++ mem")
    println("-"^80)

    for r in results
        @printf("%-8d | %-16d | %-16s | %-16d | %-16s\n",
                r.n,
                r.julia.allocs, Base.format_bytes(r.julia.memory),
                r.cpp.allocs, Base.format_bytes(r.cpp.memory))
    end

    return results
end

# Main
if abspath(PROGRAM_FILE) == @__FILE__
    use_real = "--real" in ARGS
    sizes_args = filter(x -> x != "--real", ARGS)
    sizes = isempty(sizes_args) ? [4, 8, 16, 32, 64, 128] : parse.(Int, sizes_args)

    if use_real
        run_real_benchmarks(sizes)
    else
        run_complex_benchmarks(sizes)
    end
end
