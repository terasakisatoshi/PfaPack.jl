# Native DSKR2 cancellation boundary

Independent Linux native acquisition: `/tmp/mvmc-issue176-initializer.0RKnLG`.
Command: `OPENBLAS_NUM_THREADS=1 ./probe --rank2-boundary`.
The command above describes the historical acquisition. A checkout-only
optional reproduction FROM THE AtelierArith/mvmc-rs REPOSITORY ROOT
(not from this Julia repository) is
`bash c_toolbox/issue176_rank2/reproduce.sh`; see
that toolbox README for extraction boundaries, licenses, ABI and flags.
Its smaller adapter retains the same native kernels/input without requiring
the excluded initializer toolbox. It does not replace the original provenance.
The frozen environment and complete source hashes accompany `native.txt`.
Original `fortran/dskr2.f` SHA-256:
`57bf85cfabc22055af6c95d5765200ef5c27c9b5dff074da264fc76bbf3a56f2`.
Its lines 158–176 evaluate `(A + X*TEMP1) - Y*TEMP2`.

The independent input uses upper-triangle entries
`[2^54, 1, -2^54, 2^55, 1, 2^55]`. Native factor INFO is zero
and the Pfaffian is -1. Grouping the two update products first instead
reports INFO=1: the intended nonzero final pivot is lost. The Rust regression
therefore checks factor success as an algorithm contract, not computed-bitwise
parity. Its scalar Pfaffian check allows four machine epsilons, accounting for
the two dyadic diagonal factors and extraction sign; it is not a general inverse
or factorization error budget.

Native binary SHA-256:
`0b9fb96506341c1ee247fd1ea39699516f1680c52ed248194eeec3d275618716`.
This is a standalone kernel acquisition, not a full initializer/MPI proof.
Cargo reads checked-in expectations only and does not invoke the toolbox.

Julia handoff origin: unchanged native/environment/source-sha256 files from
`tests/fixtures/issue176_initializer/rank2-boundary` in mvmc-rs. Retained
optional sources are available at
`/home/terasaki/work/atelierarith/mvmc-rs/c_toolbox/issue176_rank2/`.
The existing acquisition hashes, not an inferred current-main fixture history,
identify the native producer. No Julia-generated expectation replacement.

## Julia dependency repair validation (separate from native acquisition)

Julia1.13.1, reference Manifest09eb, numerical threads1. Actual Real RED
loaded the fresh candidate package and rebound ONLY the original2a4120d
`julia_dsktf2!` definition in that isolated process:3PASS/4FAIL, reached
INFO1 versus nativeINFO0. This is not a full oldcheckout/model proof.
Corrected fresh package GREEN:7/7PASS (direct and Float64 delegating route).
Retained receipt: `/tmp/mvmc-issue176-julia-real-validation-proof.p4r6uO`.

That acquisition later terminated83 at an ownership check: the inner timeout
created another PGID. Adapter cleanup0 covered only its leader group; remaining
owned timeout/Julia were identity-checked and explicitly TERM-cleaned, then
accessible SID/PGID scan found no remaining processes. This gap is retained
in `remaining-owned-before.json`/`remaining-owned-after.json`, not hidden.

Approved continuation used foreground timeout, without repeating RED/GREEN:
`/tmp/mvmc-issue176-julia-real-continuation-proof.CuMVwk`, terminal0.
ExistingPfaPack5699PASS plusComplex5 andReal7; outer ordinaryReal13PASS.
Source/runtime/launcher posts and owned-group-empty0. Actual loaded providers
bound to frozen files. No Fortran tests skipped; native prerequisites retained.
This validates the narrow factor/caller repair, not FSZ/SR/PhysCal/CG/MPI/Mac
closure, and does not supply a tolerance budget for downstream col13.
