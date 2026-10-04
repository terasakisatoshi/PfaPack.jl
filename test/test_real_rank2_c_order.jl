using Test, PfaPack

@testset "Real DSKR2 independent native cancellation contract" begin
    line = first(readlines(joinpath(@__DIR__, "fixtures", "issue176_real_rank2", "native.txt")))
    @test line == "rank2-boundary status=0 factor=0 inverse_calls=1 pf=-1"
    original = zeros(4,4)
    for (i,j,value) in ((1,2,2.0^54), (1,3,1.0), (1,4,-2.0^54),
                         (2,3,2.0^55), (2,4,1.0), (3,4,2.0^55))
        original[i,j] = value
        original[j,i] = -value
    end
    # Float64 Real dispatch delegates to julia_dsktf2!, not fallback BLAS ger.
    for factor! in (julia_dsktf2!, julia_zsktf2!)
        A = copy(original)
        piv = zeros(Int,4)
        @test factor!(A,piv) == 0
        @test all(p -> 1 <= p <= 4,piv) # native fixture did not record pivots
        pf = PfaPack.utu2pfa(4,A,4,piv)
        @test abs(pf - (-1.0)) <= 4eps(Float64)
    end
    # Native inverse values are diagnostics, not an inverse forward gate.
end
