using Test, PfaPack

@testset "Complex rank2 native C cancellation contract" begin
    # Independent GNU13 native ZSKTRF/child fixture, not Julia-generated values.
    # See fixtures/issue176_rank2_complex/README.md for acquisition provenance.
    root = joinpath(@__DIR__, "fixtures", "issue176_rank2_complex")
    @test readlines(joinpath(root, "native.txt")) == [
        "complex-rank2-boundary status=0 factor=0 pfa_calls=1 inverse_calls=1 pf_re=-1 pf_im=0",
        "pivots 1 2 3 4",
    ]
    A = zeros(ComplexF64, 4, 4)
    for (i, j, value) in ((1,2,2.0^54), (1,3,1.0), (1,4,-2.0^54),
                         (2,3,2.0^55), (2,4,1.0), (3,4,2.0^55))
        A[i,j] = value
        A[j,i] = -value
    end
    piv = zeros(Int, 4)
    @test julia_zsktf2!(A, piv) == 0
    @test piv == [1,2,3,4]
    pf = PfaPack.utu2pfa(4, A, 4, piv)
    # Native PF=-1; dyadic products isolate association, not a global bits gate.
    @test abs(real(pf) - (-1.0)) <= 4eps(Float64)
    @test abs(imag(pf)) <= 4eps(Float64)
end
