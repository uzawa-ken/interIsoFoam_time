# 圧力ポアソン方程式のGNN出力機能

## 概要

本ソルバー（interIsoFoam_time）は、圧力ポアソン方程式の前処理後の係数行列、右辺ベクトル、解ベクトル、残差をファイル出力する機能を持っています。これらのデータはGraph Neural Network (GNN) の学習などに活用できます。

## 出力されるデータ

### 1. 前処理後の係数行列 A（CSR形式）
**ファイル名**: `A_csr_<time>_rank<rank>.dat`

Jacobi（対角スケーリング）前処理が適用された係数行列をCSR (Compressed Sparse Row) 形式で出力します。

**前処理の定義**:
```
A_precond[i][j] = A[i][j] / diag[i]
```

これにより、対角成分は全て1.0になります。

**ファイルフォーマット**:
```
# Preconditioned matrix (Jacobi/diagonal scaling)
# A_precond[i][j] = A[i][j] / diag[i]
nRows <行数>
nCols <列数>
nnz <非ゼロ要素数>
ROW_PTR
<行ポインタ配列（スペース区切り）>
COL_IND
<列インデックス配列（スペース区切り）>
VALUES
<前処理後の値配列（スペース区切り）>
```

### 2. 前処理後の右辺ベクトル b
**ファイル名**: `b_<time>_rank<rank>.dat`

**前処理の定義**:
```
b_precond[i] = b[i] / diag[i]
```

**ファイルフォーマット**:
```
# Preconditioned RHS vector (Jacobi/diagonal scaling)
# b_precond[i] = b[i] / diag[i]
nCells <セル数>
<cellIndex> <b_precond[i]>
...
```

### 3. 解ベクトル x（時刻 n+1 の圧力）
**ファイル名**: `x_<time>_rank<rank>.dat`

圧力ポアソン方程式を解いた後の圧力場 p_rgh（時刻 n+1）を出力します。

**ファイルフォーマット**:
```
# Solution vector x (pressure p_rgh at time n+1)
nCells <セル数>
<cellIndex> <x[i]>
...
```

### 4. 残差 r
**ファイル名**: `r_<time>_rank<rank>.dat`

残差 r = Ax - b と前処理後の残差の両方を出力します。

**前処理後の残差の定義**:
```
r_precond[i] = r[i] / diag[i]
```

**ファイルフォーマット**:
```
# Residual r = Ax - b
# Columns: cellIndex  r[i]  r_precond[i]
# r_precond[i] = r[i] / diag[i] (Jacobi preconditioned)
nCells <セル数>
<cellIndex> <r[i]> <r_precond[i]>
...
```

## 有効化方法

`system/fvSolution` ファイルに以下の設定を追加します：

```
gnnOutput
{
    writePressureSystem  true;   // GNN出力を有効化
    writeInterval        1;      // 出力間隔（タイムステップ数）
}
```

## 出力ディレクトリ

出力ファイルは `<case>/gnn/` ディレクトリに保存されます：

```
<case>/
└── gnn/
    ├── A_csr_0.001_rank0.dat      # 前処理後の係数行列
    ├── b_0.001_rank0.dat          # 前処理後の右辺ベクトル
    ├── x_0.001_rank0.dat          # 解ベクトル（圧力）
    ├── r_0.001_rank0.dat          # 残差
    ├── pEqn_0.001_rank0.dat       # ノード・エッジ情報
    └── divPhi_0.001_rank0.dat     # フラックス発散
```

## 並列計算への対応

MPI並列計算時は、各プロセスが独立してファイルを出力します。ファイル名には `rank<rank番号>` が付加されます。

## 数学的背景

### 圧力ポアソン方程式

二相流ソルバーでは、PIMPLE法に基づいて以下の圧力ポアソン方程式を解きます：

```
∇·(rAUf·∇p_rgh) = ∇·(phiHbyA)
```

ここで：
- `rAUf`: 運動量方程式の逆行列係数（面補間値）
- `p_rgh`: 修正圧力（p - ρgh）
- `phiHbyA`: 修正ボリュームフラックス

### Jacobi前処理

対角スケーリング（Jacobi）前処理は最も基本的な前処理手法です：

**元の系**: Ax = b

**前処理後の系**: D⁻¹Ax = D⁻¹b

ここで D は A の対角成分からなる対角行列です。

前処理後の行列の特徴：
- 対角成分は全て 1.0
- スケーリングにより数値的安定性が向上
- 反復法の収束が改善される場合がある

## 追加の出力データ

### pEqn_*.dat（ノード・エッジ情報）

各セルとエッジの詳細情報を出力します：

**セル情報**:
- セル中心座標 (Cx, Cy, Cz)
- 対角成分 diag[i]
- 右辺 b[i]（前処理前）
- メッシュ品質メトリクス（skewness, nonOrthogonality, aspectRatio）
- 対角コントラスト、セル体積、セルサイズ、サイズジャンプ
- セルCourant数

**エッジ情報**:
- 面インデックス、所有セル、隣接セル
- 下三角成分 lower[f]、上三角成分 upper[f]

### divPhi_*.dat（フラックス発散）

フラックスの発散値を出力します：
- div(phiHbyA): 修正前フラックスの発散
- div(phi): 修正後フラックスの発散
- セル体積

## 使用例

### Python での読み込み例

```python
import numpy as np
from scipy.sparse import csr_matrix

def load_preconditioned_matrix(filename):
    """前処理後の係数行列をCSR形式で読み込む"""
    with open(filename, 'r') as f:
        lines = f.readlines()

    # ヘッダーをスキップしてサイズ情報を取得
    nRows = int(lines[2].split()[1])
    nCols = int(lines[3].split()[1])
    nnz = int(lines[4].split()[1])

    # ROW_PTR, COL_IND, VALUES を解析
    row_ptr = np.array(lines[6].split(), dtype=int)
    col_ind = np.array(lines[8].split(), dtype=int)
    values = np.array(lines[10].split(), dtype=float)

    return csr_matrix((values, col_ind, row_ptr), shape=(nRows, nCols))

def load_vector(filename):
    """ベクトルデータを読み込む"""
    data = np.loadtxt(filename, skiprows=3)
    return data[:, 1]  # 2列目が値

# 使用例
A = load_preconditioned_matrix('gnn/A_csr_0.001_rank0.dat')
b = load_vector('gnn/b_0.001_rank0.dat')
x = load_vector('gnn/x_0.001_rank0.dat')
r = np.loadtxt('gnn/r_0.001_rank0.dat', skiprows=4)
r_original = r[:, 1]
r_precond = r[:, 2]
```

## 注意事項

1. 出力ファイルのサイズは大規模メッシュでは非常に大きくなる可能性があります
2. 出力間隔 `writeInterval` を適切に設定してディスク使用量を管理してください
3. 前処理後の行列は対角成分が 1.0 であることを確認できます

## バージョン履歴

- 2024-01: 前処理後の行列・ベクトル出力機能を追加
