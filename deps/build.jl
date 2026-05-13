using Libdl

# Get the directory where this build script is located
const depsdir = @__DIR__

# Get the appropriate library extension for the current platform
libext = dlext

# Build zsktf2 library (Fortran - complex version)
libname_zsktf2 = joinpath(depsdir, "libzsktf2.$libext")
zskr2_f = joinpath(depsdir, "zskr2.f")
zsktf2_f = joinpath(depsdir, "zsktf2.f")
if Sys.isapple()
    run(`gfortran -shared -fPIC -O3 -o $libname_zsktf2 $zskr2_f $zsktf2_f -framework Accelerate`)
else
    run(`gfortran -shared -fPIC -O3 -o $libname_zsktf2 $zskr2_f $zsktf2_f -lblas -llapack`)
end

# Build dsktf2 library (Fortran - real version)
libname_dsktf2 = joinpath(depsdir, "libdsktf2.$libext")
dskr2_f = joinpath(depsdir, "dskr2.f")
dsktf2_f = joinpath(depsdir, "dsktf2.f")
if Sys.isapple()
    run(`gfortran -shared -fPIC -O3 -o $libname_dsktf2 $dskr2_f $dsktf2_f -framework Accelerate`)
else
    run(`gfortran -shared -fPIC -O3 -o $libname_dsktf2 $dskr2_f $dsktf2_f -lblas -llapack`)
end

# Build ltl2inv library (C++ — xrq-phys/Pfaffine ports, MPL-2.0)
# Note: blalink_gemmt.c is intentionally excluded; it requires the BLIS
# GEMMT extension and is not needed for the ltl2inv entry point.
# `ilaenv_lauum.cc` uses the SLOPPY_ILAENV fast-path (hard-coded
# blocksize = 64), so no separate Fortran helper is needed at link time.
libname_ltl2inv = joinpath(depsdir, "libltl2inv.$libext")
ltl2inv_cc = joinpath(depsdir, "ltl2inv.cc")
ilaenv_lauum_cc = joinpath(depsdir, "ilaenv_lauum.cc")

if Sys.isapple()
    run(`g++ -shared -fPIC -std=c++11 -O3 -D_CC_IMPL -DBLAS_EXTERNAL -I$depsdir \
        -o $libname_ltl2inv \
        $ltl2inv_cc $ilaenv_lauum_cc \
        -framework Accelerate`)
else
    run(`g++ -shared -fPIC -std=c++11 -O3 -D_CC_IMPL -DBLAS_EXTERNAL -I$depsdir \
        -o $libname_ltl2inv \
        $ltl2inv_cc $ilaenv_lauum_cc \
        -lblas -llapack`)
end

println("Built: $libname_zsktf2")
println("Built: $libname_dsktf2")
println("Built: $libname_ltl2inv")
