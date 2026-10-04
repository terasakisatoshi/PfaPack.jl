# Native complex cancellation boundary

Acquired once in Linux x86_64 container 73c57e563c61, stage
`/tmp/mvmc-issue176-rank2-complex.IbcE7c`, handle 73958 terminal exit 0.
GNU 13.3, `-O0 -ffp-contract=off`, LP64 OpenBLAS 0.3.26 Haswell,
actual numerical threads 1. Serial child check, not MPI or sampling proof.

Origin: AtelierArith/mvmc-rs, not a program shipped in this Julia repository.
Optional reproduction from that Rust repository's root:
`bash c_toolbox/issue176_rank2_complex/reproduce.sh`.
Retained source is available at
`/home/terasaki/work/atelierarith/mvmc-rs/c_toolbox/issue176_rank2_complex/`;
source revision for this handoff is main12c5bd83967c12b8d246fa81a102ada360bc19e1
plus the previously reviewed, uncommitted issue176 toolbox sources, whose exact
acquisition probe/script hashes are pinned below (not claimed present in main).
Original GPL matrix.c lines 1–21,125–179 and unchanged whole complex
Fortran kernels; source/extraction/ABI details are documented in that toolbox.
All compiler inputs were checked before compilation and after execution;
six resolved linked libraries were checked afterwards (bounded ldd closure).

Acquisition probe SHA256: 33da301fe2a6bab03cac5572dd0ffb5fc6b4af1596525947fab749f0f704afe9.
Script SHA256: 9618481eabc4536da09bf2c2213f90fe0967abeb79f7c63255499375ebcfc2c6.
Native binary SHA256: db8eceabb32cccb00474a1b6805bcd404d5d496e3b62683508c74447c0f1e069.
Native output SHA256: 7330cdf2f871eee2eccb2244a3dd683d4fba2b49db63528430cc10dbc71e2bf2.

The fixture is the unchanged native stdout, not Rust-generated expectations.
Tests require exact factor status/pivots and a four-EPSILON scaled scalar PF
comparison; no computed-bitwise acceptance or inverse forward bound.

Historical Rust test-only RED on isolated main `19a3dd29` (NOT Julia RED): run
`601b0c45-efe8-4e95-8838-6b7a8f04d7f3`, exit 100, both ordinary and
C-compatible kernels returned INFO=1 versus native INFO=0. No quotient
closure overlay was needed to expose this zero-imaginary association defect.
Retained log: `/tmp/mvmc-issue176-complex-red.log`.
