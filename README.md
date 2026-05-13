# PfaPack.jl

A Julia library for the Pfaffian and LTL decomposition of dense
skew-symmetric matrices, combining (i) pure-Julia implementations
adapted from upstream [Wimmer PfaPack](http://www.michaelwimmer.org/downloads.html)
with (ii) a C++ reference path (`xrq-phys/Pfaffine` ports) accessible
via thin ccall wrappers. Both implementation paths are unit-tested
against each other to bit-for-bit / 1e-12 tolerance.

The package provides:

- **Pfaffian computation** for skew-symmetric matrices (Parlett–Reid / LTL-based).
- **LTL decomposition** of skew-symmetric matrices (`julia_zsktf2!`, `julia_dsktf2!`, plus a `LoopVectorization`-optimised `julia_zsktf2_turbo!`).
- **Pfaffian / inverse extraction from an existing LTL form** (`utu2pfa`, `utu2inv!`, plus a C++ reference `cimpl_utu2inv!`).
- **Fortran wrappers** for the bundled PfaPack kernels (`fimpl_zsktf2!`, `fimpl_dsktf2!`).

## Status

The pure-Julia routines reproduce the bundled Fortran and C++ references to ~1e-12 across the test matrix sizes (see `test/runtests.jl`). The public API is expected to remain stable but is not yet registered in the Julia General registry.

## Native build prerequisites

`deps/build.jl` compiles three shared libraries from the bundled sources:

- `libzsktf2`, `libdsktf2` — Fortran kernels for complex/real LTL decomposition (built from `zsktf2.f`, `zskr2.f`, `dsktf2.f`, `dskr2.f`, taken byte-identically from upstream Wimmer PfaPack 2014-09).
- `libltl2inv` — C++ kernels for the `utu2inv` / `utu2pfa` family, from `xrq-phys/Pfaffine` (MPL-2.0). The build uses the SLOPPY_ILAENV fast path so no separate Fortran helper is needed.

You need:

- A Fortran compiler (`gfortran`)
- A C++11 compiler (`g++` / `clang++`)
- BLAS / LAPACK (on macOS, Apple's Accelerate framework; on Linux, `libblas`/`liblapack` or equivalent)

On Ubuntu / Debian:

```bash
sudo apt-get install -y gfortran g++ libopenblas-dev liblapack-dev
```

On macOS, the Xcode Command Line Tools plus a Homebrew `gcc` (for `gfortran`) are sufficient:

```bash
xcode-select --install
brew install gcc
```

## Install

```bash
git clone https://github.com/tmisawa/PfaPack.jl
cd PfaPack.jl
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.build()'
```

Or, from another Julia project:

```julia
using Pkg
Pkg.develop(path="/path/to/PfaPack.jl")
Pkg.build("PfaPack")
```

> `Pkg.build()` is **required** before any `using PfaPack` / `Pkg.test()`; it
> compiles `libzsktf2`, `libdsktf2`, and `libltl2inv` into `deps/`. Without
> them, the Fortran/C++ entry points (`fimpl_zsktf2!`, `fimpl_dsktf2!`,
> `cimpl_utu2inv!`) fail to load.

## Usage

### Pfaffian of a skew-symmetric matrix

```julia
using PfaPack
using LinearAlgebra

n = 4
A = zeros(n, n)
A[1, 2] = 1.0; A[1, 3] = 2.0; A[1, 4] = 3.0
A[2, 3] = 4.0; A[2, 4] = 5.0; A[3, 4] = 6.0
for i in 1:n, j in (i+1):n
    A[j, i] = -A[i, j]
end

pf = pfaffian_ltl!(copy(A))          # Parlett-Reid / LTL
```

### LTL decomposition

```julia
# Complex skew-symmetric matrix
A_c = ComplexF64.(A)
iPiv = zeros(Int, n)
info = julia_zsktf2!(A_c, iPiv)

# Real skew-symmetric matrix
A_r = copy(A)
iPiv_r = zeros(Int, n)
info_r = julia_dsktf2!(A_r, iPiv_r)
```

## Exported functions

### Pfaffian / inversion

| Function | Notes |
|---|---|
| `pfaffian_ltl!(A; overwrite_a=true)` | Pfaffian via LTL (Parlett–Reid). |
| `utu2pfa(n, A, ldA, iPiv)` | Pfaffian from an existing LTL form (Julia). |
| `utu2inv!(n, A, ldA, iPiv, vT, M, ldM)` | Inverse from LTL form (Julia). |
| `cimpl_utu2inv!(...)` | C++ reference inverse via `libltl2inv` ccall. |

### LTL decomposition

| Function | Notes |
|---|---|
| `julia_zsktf2!`, `julia_dsktf2!` | Pure-Julia LTL for complex / real matrices (upper-triangular update). |
| `julia_zsktf2_turbo!` | `LoopVectorization`+`StructArrays` SIMD variant for ComplexF64. |
| `fimpl_zsktf2!`, `fimpl_dsktf2!` | Fortran reference (bundled PfaPack kernels). |

## Test

The full test suite exercises the bundled Fortran/C++ kernels, so
`Pkg.build()` must succeed first:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.build()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Notes on numerics

- LTL decomposition: $P^T A P = L T L^T$ with $L$ unit lower-triangular and $T$ skew-tridiagonal. The Pfaffian and inverse follow from the tridiagonal factor; this avoids forming the dense LU of $A$ and preserves the skew-symmetric structure to round-off.
- All routines act on column-major Julia arrays and update only the upper-triangular part (matching the bundled Fortran kernels' `UPLO='U'` convention).

## License

This package uses a mixed-license model:

- **Primary license**: **BSD 3-Clause** — for the Julia wrapper code in `src/{PfaPack.jl, pfaffian.jl, ltl_decomposition.jl, c_wrapper.jl, fortran_wrapper.jl}` and the build/glue code `deps/build.jl`. See [`LICENSE`](LICENSE).
- **MPL-2.0 sub-components**: `src/utu2.jl` (Julia ports of `utu2pfa`, `utu2inv!`, internal `sktdsmx!`) and the C++ wrapper sources under `deps/` (`blalink.hh`, `blalink_fort.h`, `blalink_gemmt.{c,hh}`, `colmaj.hh`, `pfaffian.tcc`, `invert.tcc`, `trmmt.tcc`, `ilaenv_lauum.{cc,hh}`, `ltl2inv.cc`) — derived from [xrq-phys/Pfaffine](https://github.com/xrq-phys/Pfaffine) under the Mozilla Public License 2.0. Each file carries the standard MPL 2.0 per-file notice; the full license text is in [`LICENSE.MPL-2.0`](LICENSE.MPL-2.0).
- **Bundled upstream Wimmer Fortran kernels** (`deps/zsktf2.f`, `zskr2.f`, `dsktf2.f`, `dskr2.f`): byte-identical to upstream Wimmer PfaPack 2014-09 under the LAPACK-style modified BSD-3 (`LapackLicence`).
- **Bundled BLIS header** (`deps/blis.h`): BSD-3-Clause (BLIS upstream).
- **`pfaffian_ltl!` in `src/pfaffian.jl`**: line-by-line port of the same Wimmer release's `python/pfaffian.py`, distributed under the BSD-3-Clause wrapper license, with the upstream `LapackLicence` as the source-of-derivation license (BSD-3-compatible).

MPL 2.0 is file-scoped weak copyleft and is compatible with both the BSD-3 primary license and downstream GPL consumers (e.g. Julia-mVMC); the per-file MPL notices stay in place on the relevant files in this repository.

See [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) for the full per-file copyright/attribution mapping, upstream sources, and reproduced license texts.

## Citation

If you use this package in academic work, please cite Wimmer's PfaPack paper for the underlying Pfaffian algorithms:

> Michael Wimmer, *Algorithm 923: Efficient Numerical Computation of the Pfaffian for Dense and Banded Skew-Symmetric Matrices*, ACM Transactions on Mathematical Software **38** (2012), 30:1–30:17. <https://doi.org/10.1145/2331130.2331138>

For the C++ `ltl2inv` family (`utu2pfa` / `utu2inv` / `sktdsmx` and the
`blalink` BLIS wrappers), please also acknowledge the upstream
implementation by RuQing Xu in [xrq-phys/Pfaffine](https://github.com/xrq-phys/Pfaffine).
