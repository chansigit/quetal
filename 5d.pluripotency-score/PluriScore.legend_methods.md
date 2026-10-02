# Pluripotency gene-set score figures – legends and methods

File: `../figures/Figure1.PluriScore.MDEmap.pdf` · code: `pluri_score_mde.py`

## Figure 1G (right) – Figure legend

**(G, right) Pluripotency gene-set score on the MDE embedding.** Cells of each sample are coloured by a pluripotency score (*NANOG, FOXD3, GDF3, UTF1, TERT, PRDM14, DPPA3, DPPA5*; `scanpy.tl.score_genes`, see Methods), from 0.05 (grey) to the 99.5th percentile of all cells (0.29, blue). Dashed line, outline of all cells in the dataset. hPSC, pluripotency medium; 12h/24h Ectoderm, 12 h or 24 h of ectoderm induction; → 1d / 4d Plur. media, induced cells returned to pluripotency medium for 1 or 4 days. Arrow, 12 h-induced cells that re-occupy the hPSC region after 4 days in pluripotency medium.

## Methods

**Embedding.** All panels use the canonical MDE embedding of the 27-sample H1 dataset (`adata_merged.250505-canonical.h5ad`). Briefly, after QC and normalisation, an scVI latent space was learned with nuisance covariates (stress, mitochondrial, sex and cell-cycle programmes) regressed out, batch-corrected with Harmony (pool and genetic background), and embedded in two dimensions with minimum-distortion embedding (MDE, pymde) anchored on the twelve control samples; perturbed samples were projected onto this anchored embedding. The dashed outline in every panel is a concave hull (alpha shape, α = 3.5, 0.03 buffer; boundary resampled to 1,500 points and Gaussian-smoothed, σ = 18 points) of 8,000 randomly chosen cells from all 27 samples, and is identical in all panels so that per-sample occupancy can be compared.

**Pluripotency score.** Gene-set scores were computed with `scanpy.tl.score_genes` (scanpy 1.12.4) on the log1p library-size-normalised expression matrix, jointly for all 123,822 cells of the 27 samples, with the gene list *NANOG, FOXD3, GDF3, UTF1, TERT, PRDM14, DPPA3, DPPA5*, `ctrl_size = 50`, `n_bins = 25`, `random_state = 0`. Following the Seurat `AddModuleScore` procedure, all genes are binned into 25 bins by mean expression, 50 control genes are sampled from the bin of each gene in the set, and the score of a cell is the mean expression of the gene set minus the mean expression of the control genes. The score is therefore in units of log1p-normalised expression, is centred at 0 for cells in which the gene set is expressed no higher than expression-matched background genes, and is positive when the pluripotency genes are coordinately up-regulated. Scoring all samples together places Figure 1 and Figure 2 panels on one scale.

**Visualisation.** For each sample, cells were plotted at their MDE coordinates (point size 1.2 pt, rasterised) and coloured with a three-colour linear map (#d1cfd4 → #9577e5 → #1206f5) spanning a score of 0.05 (lower bound, chosen so that cells with no or marginal enrichment over background share the base grey) to the 99.5th percentile of the score over all cells (0.29); scores below or above these bounds are clipped to the end colours. Cells were drawn in ascending order of score. The per-cell scores and MDE coordinates used for the figure are stored in `pluri_score.cache.parquet`.

---

# Figure 2B – Control vs *JARID2*-CRISPRi

Files: `../figures/Figure2.PluriScore.MDEmap.pdf` (top, code `pluri_score_mde_fig2.py`), `../figures/Figure2.PluriScore.violin.pdf` + `Figure2.PluriScore.violin.stats.csv` (this directory) (bottom, code `pluri_score_violin_fig2.py`)

## Figure legend

**(B) *JARID2*-deficient ectoderm cells reacquire a pluripotency-like expression programme.** Top, MDE embedding of control (left) and *JARID2*-CRISPRi (right) cells (hPSC, ectoderm, and 24 h ectoderm returned to pluripotency medium for 1 or 4 days), coloured by the pluripotency gene-set score (same score as in G; colour range 0.05 to the 99.5th percentile of the cells shown, 0.34). Arrows, inferred trajectories. Bottom, distribution of the pluripotency score per condition (violins truncated at the data range, equal width). White dots, medians; Δ, difference of medians (*JARID2*-CRISPRi − control). ***, P < 0.001, two-sided Mann–Whitney U test within each condition.

## Methods (additions)

**Figure 2 panels.** The MDE coordinates, score and colour map are those described above; the two MDE panels share axis limits, and the colour range was set to 0.05–0.34 (99.5th percentile of the score over the eight Figure 2 samples). Violins were drawn with seaborn (`cut = 0`, `density_norm = 'width'`, no inner marks) and coloured as the sample map in Figure 2. For each condition, per-cell scores of *JARID2*-CRISPRi and control cells were compared with a two-sided Mann–Whitney U test, and the difference of medians (Δ) is reported as effect size (n = 1,224–11,926 cells per group; all P < 1e-30, see `Figure2.PluriScore.violin.stats.csv`). Because these are per-cell tests on thousands of cells, effect sizes (median difference 0.03–0.09 score units) are more informative than the P values.

---

# Figure 2 (supplementary grid) – per-sample score maps, Control vs *JARID2*-CRISPRi

File: `../figures/Figure2.PluriScore.MDEmap.grid.pdf` · code `pluri_score_mde_fig2_grid.py`

## Figure legend

**Per-sample pluripotency score on the MDE embedding, control (top) and *JARID2*-CRISPRi (bottom).** Columns: hPSC, 24 h ectoderm, and 24 h ectoderm returned to pluripotency medium for 1 or 4 days. Colour, pluripotency score from 0.10 (grey) to 0.30 (blue). Dashed line, outline of all cells in the dataset. Numbers, percentage of cells of the sample with a score above 0.15.

## Methods (addition)

The grid uses the same embedding, hull and score as Figure 1G. The colour range was fixed at 0.10–0.30 so that cells re-acquiring a high score are resolved against the bulk of low-scoring cells; the percentage of cells with score > 0.15 is given per panel as a threshold-based summary (full threshold table in `pluri_score.threshold_fractions.csv`). In the 4-day condition this is 5.0 % (control) vs 16.6 % (*JARID2*-CRISPRi; 3.3-fold), and 27 % vs 52 % of cells lie within the hPSC region of the embedding, whereas the median score of cells already inside the hPSC region is similar between the two (0.063 vs 0.072).
