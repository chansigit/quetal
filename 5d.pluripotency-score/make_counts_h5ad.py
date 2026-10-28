"""Rewrite the canonical object with raw counts as X (CSR, compressed) and no log-normalised matrix.

    python make_counts_h5ad.py <in.h5ad> <out.h5ad>

Keeps obs / var / obsm / obsp / uns / varm untouched. layers['counts'] becomes X (float32,
integer-valued, CSR); the old log1p X and the duplicate layer are dropped. Written with gzip
compression. Prints a verification (shapes, nnz, element-wise equality with the source counts).
"""

import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp

src, dst = Path(sys.argv[1]), Path(sys.argv[2])

a = ad.read_h5ad(src, backed='r')   # X (log matrix) stays on disk; only layers/obs/obsm are loaded
counts = a.layers["counts"]
assert sp.issparse(counts), "counts layer is not sparse"
counts = sp.csr_matrix(counts, dtype=np.float32)
assert (counts.data == np.round(counts.data)).all(), "counts are not integer-valued"

b = ad.AnnData(
    X=counts,
    obs=a.obs,
    var=a.var,
    obsm=dict(a.obsm),
    obsp=dict(a.obsp),
    varm=dict(a.varm),
    uns=dict(a.uns),
)
b.uns["X_is"] = (
    "raw UMI counts (CSR); log-normalised matrix dropped 2026-10-28, recompute with sc.pp.normalize_total + sc.pp.log1p"
)
a.file.close(); del a
b.write_h5ad(dst, compression="gzip", compression_opts=4)

# verify against the source file without loading it all again
with h5py.File(src) as f, h5py.File(dst) as g:
    for k in ["data", "indices", "indptr"]:
        s, d = f["layers/counts"][k], g["X"][k]
        assert s.shape == d.shape, k
        n = s.shape[0]
        for i in range(0, n, 50_000_000):
            assert np.array_equal(s[i : i + 50_000_000], d[i : i + 50_000_000]), (
                f"{k} mismatch at {i}"
            )
    print("verified: X of", dst.name, "== layers/counts of", src.name)
    print("obsm:", list(g["obsm"].keys()))
print(
    f"{src.name}: {src.stat().st_size / 1e9:.2f} GB  ->  {dst.name}: {dst.stat().st_size / 1e9:.2f} GB"
)
