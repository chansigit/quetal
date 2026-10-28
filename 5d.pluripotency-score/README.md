# Pluripotency gene-set score figures

Run from `workflow-250505/` with the dl2025 venv (`source ~/pp`), e.g. `python 5d.pluripotency-score/pluri_score_mde.py`. PDFs go to the shared `figures/`.

| Script | Output | Notes |
|---|---|---|
| `pluri_score_mde.py` | `figures/Figure1.PluriScore.MDEmap.pdf` | Fig 1G right part. First run reads `adata_merged.250505-canonical.h5ad`, computes `sc.tl.score_genes` on all 27 samples and writes `5d.pluripotency-score/pluri_score.cache.parquet`; later runs use the cache. Options: `--recompute`, `--order random`, `--vmin`, `--vmax-pct`, `--tag`. |
| `pluri_score_mde_fig2.py` | `figures/Figure2.PluriScore.MDEmap.pdf` | Fig 2B top: Control vs JARID2-CRISPRi MDE panels coloured by score. Needs the cache. |
| `pluri_score_mde_fig2_grid.py` | `figures/Figure2.PluriScore.MDEmap.grid.pdf` | Fig 2 supplement: 2×4 grid, Control (top) vs JARID2-CRISPRi (bottom), Fig 1G-style hull panels; colour range 0.10–0.30, per-panel % of cells with score > 0.15 (`--vmin --vmax --thr --decimals`). `pluri_score.threshold_fractions.csv` = fractions above a series of thresholds per sample. Needs the cache. |
| `pluri_score_violin_fig2.py` | `figures/Figure2.PluriScore.violin.pdf`, `.stats.csv` | Fig 2B bottom: violins, medians, Δ median, two-sided Mann-Whitney U. Needs the cache. |

Gene set: NANOG, FOXD3, GDF3, UTF1, TERT, PRDM14, DPPA3, DPPA5. Legends and methods: `PluriScore.legend_methods.md` (this directory).

## Compact counts object

`make_counts_h5ad.py <in> <out>` rewrites the canonical object with raw UMI counts as `X` (CSR, gzip) and drops the log-normalised matrix and the duplicate `layers['counts']`; obs/var/obsm/obsp/uns are kept. 6.56 GB -> 1.02 GB, verified element-wise against the source. `check_counts_h5ad.py <h5ad>` re-opens the result and runs sanity checks. Both need ~8-12 GB RAM: run on a compute node (`srun -p dev --mem=16G`), not a login node. Output lives in `/oak/stanford/projects/eca/sc-data/yimiao-quetal/adata_merged.250505-canonical.counts.h5ad`. Note: the original log-normalised X cannot be reproduced exactly from these counts (normalisation was done before nuisance genes were removed from `var`).
