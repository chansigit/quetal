import h5py, anndata as ad, numpy as np, scipy.sparse as sp, sys
p=sys.argv[1]
with h5py.File(p,'r') as f:
    print('HDF5 opens; top-level keys:', list(f.keys()))
    print('X encoding:', dict(f['X'].attrs), '| compression:', f['X/data'].compression)
a=ad.read_h5ad(p)
print('AnnData:', a.shape, '| X type:', type(a.X).__name__, a.X.dtype, '| nnz:', a.X.nnz)
assert sp.isspmatrix_csr(a.X)
print('integer-valued:', bool((a.X.data==np.round(a.X.data)).all()), '| min/max:', a.X.data.min(), a.X.data.max())
print('obs cols:', len(a.obs.columns), '| samples:', a.obs["sample"].nunique(), '| var:', a.n_vars, '| NANOG in var:', 'NANOG' in a.var_names)
print('obsm:', sorted(a.obsm.keys()))
print('X_mde shape:', a.obsm['X_mde'].shape, '| uns X_is:', str(a.uns['X_is'])[:60])
tot=np.asarray(a.X.sum(1)).ravel(); print('total counts per cell: median', int(np.median(tot)), 'min', int(tot.min()))
hp=a[a.obs['sample']=='Ctrl-hPSC.p1']; g=hp[:, 'NANOG'].X.toarray()
print('NANOG mean count in Ctrl-hPSC.p1:', round(float(g.mean()),3), '| frac expressing:', round(float((g>0).mean()),3))
print('VALID')
