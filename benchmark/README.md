# Benchmark

PfaPack.jl のパフォーマンス比較ベンチマーク

## ベンチマークスクリプト

### 1. LTL分解 (`benchmark_ltl.jl`)

LTL分解の実装比較

#### 比較対象

1. **`julia_dsktf2!`** - 最適化されたJulia実装（Float64用、上三角のみ更新）
2. **`julia_zsktf2!`** - 純粋なJulia実装（ComplexF64用）
3. **`julia_zsktf2_turbo!`** - LoopVectorization.jl + StructArrays.jl による最適化実装（ComplexF64用）
4. **`fimpl_dsktf2!`** - Fortran実装（Float64用）
5. **`fimpl_zsktf2!`** - Fortran実装（ComplexF64用）

#### 実行方法

```bash
# Real matrices (Julia vs Fortran comparison)
julia --project=. benchmark/benchmark_ltl.jl

# Complex matrices (Julia vs Turbo vs Fortran comparison)
julia --project=. benchmark/benchmark_ltl.jl --complex

# カスタムサイズを指定
julia --project=. benchmark/benchmark_ltl.jl 4 8 16 32 64 128 256 512 1024
julia --project=. benchmark/benchmark_ltl.jl --complex 4 8 16 32 64 128 256 512 1024
```

#### 結果例 (Apple M2 Max)

**Float64 (Julia vs Fortran):**
```
Size   | Julia (ms)   | Fortran (ms) | F/J      | Match
-------------------------------------------------------
4      | 0.000025     | 0.000060     | 2.45     | ✓
8      | 0.000135     | 0.000319     | 2.36     | ✓
16     | 0.000579     | 0.000897     | 1.55     | ✓
32     | 0.002893     | 0.003072     | 1.06     | ✓
64     | 0.015145     | 0.013792     | 0.91     | ✓
128    | 0.097604     | 0.093916     | 0.96     | ✓
256    | 0.999229     | 1.162229     | 1.16     | ✓
512    | 7.619313     | 12.330958    | 1.62     | ✓
1024   | 54.638146    | 98.874042    | 1.81     | ✓
```

**ComplexF64 (Julia vs Turbo vs Fortran):**
```
Size   | Julia (ms)   | Turbo (ms)   | Fortran (ms) | F/J      | F/T      | T/J
-------------------------------------------------------------------------------------
4      | 0.000038     | 0.000107     | 0.000105     | 2.75     | 0.98     | 2.81
8      | 0.000215     | 0.000295     | 0.000281     | 1.31     | 0.95     | 1.37
16     | 0.001308     | 0.001408     | 0.001007     | 0.77     | 0.72     | 1.08
32     | 0.008167     | 0.006946     | 0.004774     | 0.58     | 0.69     | 0.85
64     | 0.050625     | 0.037521     | 0.031750     | 0.63     | 0.85     | 0.74
128    | 0.446854     | 0.295000     | 0.308229     | 0.69     | 1.04     | 0.66
256    | 3.828187     | 2.492396     | 3.522125     | 0.92     | 1.41     | 0.65
512    | 31.885375    | 17.527729    | 27.834854    | 0.87     | 1.59     | 0.55
1024   | 255.115062   | 128.752125   | 188.834666   | 0.74     | 1.47     | 0.50
```

**凡例:**
- **Julia (Float64)** = `julia_dsktf2!` (最適化Julia実装、上三角のみ更新)
- **Julia (ComplexF64)** = `julia_zsktf2!` (純粋Julia実装、上三角のみ更新)
- **Turbo (ComplexF64)** = `julia_zsktf2_turbo!` (LoopVectorization.jl + StructArrays.jl)
- **Fortran (Float64)** = `fimpl_dsktf2!` mode='N' (完全なtridiagonal化)
- **Fortran (ComplexF64)** = `fimpl_zsktf2!` mode='N' (完全なtridiagonal化)
- **F/J** = Fortran time / Julia time (>1 = Julia が高速)
- **F/T** = Fortran time / Turbo time (>1 = Turbo が高速)
- **T/J** = Turbo time / Julia time (<1 = Turbo が高速)

**考察:**

- **Float64**: 小さい行列（n ≤ 32）では Julia が 1.1〜2.5 倍高速。中間サイズ（n = 64〜128）では Fortran がわずかに高速（0.91〜0.96 倍）。大きい行列（n ≥ 256）では Julia が 1.2〜1.8 倍高速

- **ComplexF64 (julia_zsktf2!)**: 小さい行列（n ≤ 8）では Julia が 1.1〜2.0 倍高速。大きい行列（n ≥ 16）では Fortran が 1.1〜1.7 倍高速

- **ComplexF64 (julia_zsktf2_turbo!)**:
  - n ≥ 32 で純粋Julia実装より 15〜50% 高速
  - **n ≥ 128 で Fortran より高速！**
    - n=128: Fortran の 1.04 倍高速
    - n=256: Fortran の 1.41 倍高速
    - n=512: Fortran の 1.59 倍高速
    - n=1024: Fortran の 1.47 倍高速
  - n < 32 では StructArray へのコピーオーバーヘッドで遅くなる

**最適化のポイント:**

- 元の Julia 実装では `BLAS.ger!` を2回呼び出し、全行列 n×n を更新していた
- `julia_dsktf2!` / `julia_zsktf2!` では Fortran DSKR2/ZSKR2 と同様に上三角部分のみを更新
- Float64 では大規模行列で Julia が高速
- **`julia_zsktf2_turbo!` の最適化手法:**
  1. StructArrays.jl で ComplexF64 行列を実部・虚部の分離した連続配列に変換
  2. LoopVectorization.jl の `@turbo` マクロで内側ループを SIMD ベクトル化
  3. 複素数演算を実部・虚部の Float64 演算に分解することで効率的な SIMD 実行が可能に

---

### 2. 逆行列計算 (`benchmark_utu2inv.jl`)

`utu2inv!` (Julia) と `cimpl_utu2inv!` (C++) の比較。`utu2inv!` は `src/utu2.jl` (MPL-2.0、xrq-phys/Pfaffine の Julia port)、`cimpl_utu2inv!` は `libltl2inv.dylib` (同じく xrq-phys/Pfaffine の C++ 実装) への ccall。

#### 比較対象

1. **`utu2inv!`** - 純粋なJulia実装
2. **`cimpl_utu2inv!`** - C++実装（ccall経由）

#### 実行方法

```bash
# Complex matrices (default)
julia --project=. benchmark/benchmark_utu2inv.jl

# Real matrices
julia --project=. benchmark/benchmark_utu2inv.jl --real

# カスタムサイズを指定
julia --project=. benchmark/benchmark_utu2inv.jl 4 8 16 32 64 128 256 512 1024
```

#### 結果例 (Apple M2 Max)

**Float64:**
```
Size     | Julia (ms)     | C++ (ms)       | C++/J      | Match
----------------------------------------------------------------------
4        | 0.000331       | 0.000428       | 1.29       | ✓
8        | 0.000552       | 0.000756       | 1.37       | ✓
16       | 0.001329       | 0.001687       | 1.27       | ✓
32       | 0.014521       | 0.008130       | 0.56       | ✓
64       | 0.038042       | 0.030083       | 0.79       | ✓
128      | 0.170833       | 0.140187       | 0.82       | ✓
256      | 0.891854       | 0.880375       | 0.99       | ✓
512      | 5.199438       | 5.574104       | 1.07       | ✓
1024     | 24.265145      | 29.434895      | 1.21       | ✓
```

**ComplexF64:**
```
Size     | Julia (ms)     | C++ (ms)       | C++/J      | Match
----------------------------------------------------------------------
4        | 0.000403       | 0.000800       | 1.98       | ✓
8        | 0.000736       | 0.001571       | 2.13       | ✓
16       | 0.002435       | 0.005437       | 2.23       | ✓
32       | 0.021729       | 0.026438       | 1.22       | ✓
64       | 0.054312       | 0.120417       | 2.22       | ✓
128      | 0.230209       | 0.553875       | 2.41       | ✓
256      | 1.617104       | 2.978958       | 1.84       | ✓
512      | 9.521250       | 17.157500      | 1.80       | ✓
1024     | 44.885688      | 104.928500     | 2.34       | ✓
```

**凡例:**
- **Julia** = `utu2inv!` (純粋Julia実装)
- **C++** = `cimpl_utu2inv!` (C++ via ccall)
- **C++/J** = C++ time / Julia time (>1 = Julia が高速)

**考察:**
- **ComplexF64**: Julia 実装が全サイズで C++ より 1.2〜2.4 倍高速
- **Float64**: 小さい行列（n ≤ 16）では Julia が 1.3〜1.4 倍高速、中間サイズ（n = 32〜256）では C++ が 1.0〜1.8 倍高速、大きい行列（n ≥ 512）では Julia が 1.1〜1.2 倍高速
- Julia 実装はほぼゼロアロケーション（n ≥ 512 でわずかなアロケーション）

## 事前準備

BenchmarkTools.jl が必要です：

```julia
using Pkg
Pkg.add("BenchmarkTools")
```

## REPL から実行

```julia
using Pkg
Pkg.activate(".")

include("benchmark/benchmark_ltl.jl")
run_real_benchmarks([4, 8, 16, 32, 64, 128, 256, 512, 1024])
run_complex_benchmarks([4, 8, 16, 32, 64, 128, 256, 512, 1024])

include("benchmark/benchmark_utu2inv.jl")
run_real_benchmarks([4, 8, 16, 32, 64, 128, 256, 512, 1024])
run_complex_benchmarks([4, 8, 16, 32, 64, 128, 256, 512, 1024])
```
