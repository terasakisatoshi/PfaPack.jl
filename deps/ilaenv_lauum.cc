/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "blalink.hh"
#include "ilaenv_lauum.hh"

/*
 * Hard-coded LAUUM blocksize. On every BLAS/LAPACK target this package
 * supports, LAPACK ILAENV(ispec=1, xLAUUM, "", ...) returns 64 for
 * x = S, D, C, Z, so returning 64 directly is equivalent to calling
 * ILAENV and avoids the Fortran-string-handling subtleties of doing so
 * from C++.
 */
#define EXPANDMAC(cctype, name)                            \
template <> int ilaenv_lauum<cctype>(uplo_t uplo, int n)   \
{                                                          \
    return 64;                                             \
}

EXPANDMAC( float,    SLAUUM )
EXPANDMAC( double,   DLAUUM )
EXPANDMAC( ccscmplx, CLAUUM )
EXPANDMAC( ccdcmplx, ZLAUUM )
#undef EXPANDMAC
