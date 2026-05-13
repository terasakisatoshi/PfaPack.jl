# Third-Party Licenses

This package is distributed under a mixed-license model:

- **Primary license** (`LICENSE`): BSD-3-Clause, applies to Julia
  wrapper code in `src/{PfaPack.jl, pfaffian.jl, ltl_decomposition.jl,
  c_wrapper.jl, fortran_wrapper.jl}` and to the build/glue code
  `deps/build.jl`.
- **MPL-2.0 sub-components**: `src/utu2.jl` plus the C++ wrapper
  sources under `deps/`, all originating from the upstream
  [xrq-phys/Pfaffine](https://github.com/xrq-phys/Pfaffine) project
  by RuQing Xu (see "xrq-phys / Pfaffine" section below).
- **Bundled Wimmer PfaPack Fortran kernels** under `deps/`: retain
  Wimmer's upstream LAPACK-style modified BSD-3 from `LapackLicence`
  (see "Wimmer PfaPack" section).
- **Bundled BLIS header**: BLIS upstream BSD-3-Clause (see "BLIS" section).
- **Derivative Julia code in `src/pfaffian.jl`** (`pfaffian_ltl!`):
  line-by-line port of upstream Wimmer Python `pfaffian.py`,
  redistributed under the package's primary BSD-3 license. The
  upstream `LapackLicence` is BSD-3-compatible.

## Per-file license mapping

| File(s) | Upstream / author | License |
|---|---|---|
| `src/PfaPack.jl`, `src/pfaffian.jl`, `src/ltl_decomposition.jl`, `src/c_wrapper.jl`, `src/fortran_wrapper.jl`, `deps/build.jl` | Satoshi Terasaki, Takahiro Misawa | BSD-3-Clause |
| `src/utu2.jl` | RuQing Xu via [xrq-phys/Pfaffine](https://github.com/xrq-phys/Pfaffine) C++ originals (Julia port) | **MPL-2.0** |
| `deps/blalink.hh`, `deps/blalink_fort.h`, `deps/colmaj.hh`, `deps/blalink_gemmt.{c,hh}` | RuQing Xu, [xrq-phys/Pfaffine](https://github.com/xrq-phys/Pfaffine) | **MPL-2.0** |
| `deps/pfaffian.tcc`, `deps/invert.tcc`, `deps/trmmt.tcc`, `deps/ilaenv_lauum.{cc,hh}`, `deps/ltl2inv.cc` | RuQing Xu, originally in his `xrq-phys/mVMC` `ltl-inv` branch (2022-03), propagated to `issp-center-dev/mVMC` `src/ltl2inv/` | **MPL-2.0** |
| `deps/dskr2.f`, `deps/zskr2.f`, `deps/dsktf2.f`, `deps/zsktf2.f` | Michael Wimmer, PfaPack 2014-09 | LAPACK-style modified BSD-3 (`LapackLicence`) |
| `deps/blis.h` | BLIS project (UT Austin, HPE, AMD) | BSD-3-Clause |

The MPL-2.0 license text is reproduced in [`LICENSE.MPL-2.0`](LICENSE.MPL-2.0).
MPL-2.0 is file-scoped weak copyleft and is compatible with both BSD-3
(primary license) and GPL (downstream consumers); MPL §3.3 explicitly
permits the file-by-file combination used here.

**Default license for unlisted files**: Any file tracked in this
repository that is not explicitly listed in the table above (for
example `Project.toml`, `.gitignore`, the test suite under `test/`,
the benchmark scripts and notes under `benchmark/`, and the top-level
`README.md`, `LICENSE`, `LICENSE.MPL-2.0`, and `THIRD_PARTY_LICENSES.md`
themselves) is distributed under the package's primary BSD-3-Clause
license.

---

## xrq-phys / Pfaffine — Mozilla Public License 2.0

- **Upstream author**: RuQing Xu ([@xrq-phys](https://github.com/xrq-phys))
- **Upstream repository (primary)**: <https://github.com/xrq-phys/Pfaffine>
- **License**: Mozilla Public License 2.0 (per-file headers, root `LICENSE`)

### Provenance of the bundled files

- `deps/blalink.hh`, `deps/blalink_fort.h`, `deps/colmaj.hh`: taken
  directly from `xrq-phys/Pfaffine` `src/`. `colmaj.hh` is the
  upstream `src/colmaj.tcc` renamed and with `dim_t` replaced by `int`.
- `deps/pfaffian.tcc`, `deps/invert.tcc`, `deps/trmmt.tcc`,
  `deps/ilaenv_lauum.{cc,hh}`, `deps/ltl2inv.cc`,
  `deps/blalink_gemmt.{c,hh}`: authored by RuQing Xu and first
  published in his `xrq-phys/mVMC` `ltl-inv` branch on 2022-03-11 /
  03-15 / 03-21, then propagated to the upstream `issp-center-dev/mVMC`
  tree under `src/ltl2inv/`. They preserve the same MPL-2.0 per-file
  headers across every downstream copy.

### Derivative Julia code in this package

`src/utu2.jl` contains direct Julia ports of the following C++ routines
from xrq-phys/Pfaffine, preserving the original MPL-2.0 license per
MPL §1.10 ("Modifications"):

| Julia function | Upstream C++ source |
|---|---|
| `utu2pfa` | `deps/pfaffian.tcc: utu2pfa()` |
| `utu2inv!` | `deps/invert.tcc:76-118` (`utu2inv` template) |
| `sktdsmx!` (internal) | `deps/invert.tcc: sktdsmx()` |

### Per-file MPL-2.0 notice

Every MPL-2.0 file in this repository carries the standard upstream
notice:

> This Source Code Form is subject to the terms of the Mozilla Public
> License, v. 2.0. If a copy of the MPL was not distributed with this
> file, You can obtain one at https://mozilla.org/MPL/2.0/.

---

## Wimmer PfaPack — modified BSD-3-Clause (LAPACK-style)

- **Bundled Fortran source**: `deps/dskr2.f`, `deps/zskr2.f`, `deps/dsktf2.f`,
  `deps/zsktf2.f`
- **Derivative Julia code**: `src/pfaffian.jl` — the `pfaffian_ltl!` function
  is a line-by-line Julia port of the Parlett–Reid LTL implementation in
  upstream `python/pfaffian.py` (function `pfaffian_LTL`, lines 247–308 of
  the 2014-09 release). Per-line cross-references to the Python source are
  retained in the Julia comments. No upstream Python source is bundled.
- **Author**: Michael Wimmer
- **Upstream release**: PfaPack 2014-09 release (downloaded from
  <http://www.michaelwimmer.org/downloads.html>)
- **Verification**: The four bundled `.f` files are byte-identical to the
  corresponding files under `pfapack/fortran/` in the upstream release.

The accompanying upstream paper is:

> Michael Wimmer, *Algorithm 923: Efficient Numerical Computation of the
> Pfaffian for Dense and Banded Skew-Symmetric Matrices*, ACM Transactions
> on Mathematical Software **38** (2012), 30:1–30:17.
> <https://doi.org/10.1145/2331130.2331138>

The upstream PfaPack package distributes the Fortran, C interface, Python,
Matlab, and Mathematica implementations together under a single
`LapackLicence` file at the package root (the only license file in the
upstream tarball). The terms reproduced below from that file therefore cover
both the bundled Fortran sources and the upstream Python source from which
`pfaffian_ltl!` was ported.

The upstream `LapackLicence` file (the sole license file at the root of the
PfaPack 2014-09 release, reproduced verbatim) states:

> The FORTRAN implementation is derived in parts on the reference
> implementation of LAPACK and BLAS on www.netlib.org.
> LAPACK is subject to the following license:
>
> ======================================================================
>
> Copyright (c) 1992-2010 The University of Tennessee and The University
>                         of Tennessee Research Foundation.  All rights
>                         reserved.
> Copyright (c) 2000-2010 The University of California Berkeley. All
>                         rights reserved.
> Copyright (c) 2006-2010 The University of Colorado Denver.  All rights
>                         reserved.
>
> $COPYRIGHT$
>
> Additional copyrights may follow
>
> $HEADER$
>
> Redistribution and use in source and binary forms, with or without
> modification, are permitted provided that the following conditions are
> met:
>
> - Redistributions of source code must retain the above copyright
>   notice, this list of conditions and the following disclaimer.
>
> - Redistributions in binary form must reproduce the above copyright
>   notice, this list of conditions and the following disclaimer listed
>   in this license in the documentation and/or other materials
>   provided with the distribution.
>
> - Neither the name of the copyright holders nor the names of its
>   contributors may be used to endorse or promote products derived from
>   this software without specific prior written permission.
>
> The copyright holders provide no reassurances that the source code
> provided does not infringe any patent, copyright, or any other
> intellectual property rights of third parties.  The copyright holders
> disclaim any liability to any recipient for claims brought against
> recipient by any third party for infringement of that parties
> intellectual property rights.
>
> THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
> "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
> LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
> A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
> OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
> SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
> LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
> DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
> THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
> (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
> OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
>
> ======================================================================

---

## BLIS — BSD-3-Clause

- **Bundled source**: `deps/blis.h` (a cut-down header excerpt from upstream
  BLIS v0.8.0, containing only the type definitions and macros required by
  the `blalink_gemmt.{c,hh}` wrappers).
- **Upstream**: <https://github.com/flame/blis>

The per-file header in `deps/blis.h` carries the upstream notice:

> Copyright (C) 2014, The University of Texas at Austin
> Copyright (C) 2016, Hewlett Packard Enterprise Development LP
> Copyright (C) 2020, Advanced Micro Devices, Inc.
>
> Redistribution and use in source and binary forms, with or without
> modification, are permitted provided that the following conditions are
> met:
>  - Redistributions of source code must retain the above copyright
>    notice, this list of conditions and the following disclaimer.
>  - Redistributions in binary form must reproduce the above copyright
>    notice, this list of conditions and the following disclaimer in the
>    documentation and/or other materials provided with the distribution.
>  - Neither the name(s) of the copyright holder(s) nor the names of its
>    contributors may be used to endorse or promote products derived
>    from this software without specific prior written permission.
>
> THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
> "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
> LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
> A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
> HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
> SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
> LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
> DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
> THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
> (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
> OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
