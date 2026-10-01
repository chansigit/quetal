"""Pluripotency gene-set score on the canonical MDE, laid out like fig-to-refine1.png.

Replaces the single-gene NANOG panels (5c.gene-scatter-mde.ipynb) with an
8-gene score (sc.tl.score_genes, Seurat AddModuleScore style). Embedding (X_mde), hull and
point style are identical to 5c. Score is computed once on all 27 canonical samples.

Figure 1 (H1 differentiation series) follows the layout/annotations of fig-to-refine1.png:
  right part only: per-sample score panels: hPSC | 12h row (Ecto, ->1d, ->4d) / 24h row (same),
  two quarter-circle arrows hPSC -> rows, colourbar, 'Reversion to hPSC-like state' callout.

First run reads the 6 GB h5ad and writes a small cache (5d.pluripotency-score/pluri_score.cache.parquet);
later runs only need the cache, so layout tweaks take seconds. The hull is recomputed from the
cached coordinates each run (all 27 samples, as in 5c).

    python 5d.pluripotency-score/pluri_score_mde.py            # uses cache if present
    python 5d.pluripotency-score/pluri_score_mde.py --recompute
    python 5d.pluripotency-score/pluri_score_mde.py --order random   # shuffled draw order instead of high-score-on-top
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans']
matplotlib.rcParams['pdf.fonttype'] = 42   # editable text in Illustrator
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.spatial import Delaunay
from shapely.geometry import MultiLineString
from shapely.ops import polygonize, unary_union

HERE = Path(__file__).resolve().parent.parent      # workflow-250505/
THIS = Path(__file__).resolve().parent             # 5d.pluripotency-score/
H5AD = HERE / 'adata_merged.250505-canonical.h5ad'
OUT = HERE / 'figures'
CACHE = THIS / 'pluri_score.cache.parquet'

PLURI_GENES = ['NANOG', 'FOXD3', 'GDF3', 'UTF1', 'TERT', 'PRDM14', 'DPPA3', 'DPPA5']
SCORE = 'pluri_score'
SCORE_LABEL = 'Pluripotency score'

# ---- Figure 1: sample -> treatment label (0.5 = 12 h, E = ectoderm, T = 1 d / 4T = 4 d back in plur. media)
FIG1_SAMPLES = ['H1-hPSC.p6', 'H1-0.5E', 'H1-0.5ET', 'H1-0.5E4T', 'H1-E', 'H1-ET', 'H1-E4T.1']
PANEL_TITLE = {
    'H1-hPSC.p6': 'hPSC',
    'H1-0.5E':    '12h Ectoderm',
    'H1-0.5ET':   '12h Ecto.→\n1d Plur. media',
    'H1-0.5E4T':  '12h Ecto.→\n4d Plur. media',
    'H1-E':       '24h Ectoderm',
    'H1-ET':      '24h Ecto.→\n1d Plur. media',
    'H1-E4T.1':   '24h Ecto.→\n4d Plur. media',
}
TITLE_COLOR = '#8b2fb0'

# ---- hull / colour parameters (hull tightened vs 5c) ----
ALPHA, BUFFER, SMOOTH_RES, GAUSS_SIGMA, SUBSAMPLE = 3.5, 0.03, 128, 18, 8000   # hull: tight + smoothed (5c: 2.0, 0.15, sigma 3)
CMAP = mcolors.LinearSegmentedColormap.from_list('custom_cmap', ['#d1cfd4', '#9577e5', '#1206f5'])
VMIN = 0.05             # score below this is drawn in the base grey (score value, not a percentile)
VMAX_PERCENTILE = 99.5  # colour saturates at this percentile of the score over all cells (--vmax-pct / --vmin override)
POINT_ORDER = 'score'   # 'score' = high scores drawn on top (as 5c); 'random' = shuffled (--order random)


# ---------------- concave hull (verbatim from 5c) ----------------
def concave_hull(points, alpha=2.0, buffer_dist=0.15, resolution=128):
    pts = np.array(points)
    tri = Delaunay(pts)
    edges = set()
    for simplex in tri.simplices:
        pa, pb, pc = pts[simplex]
        a, b, c = np.linalg.norm(pa - pb), np.linalg.norm(pb - pc), np.linalg.norm(pc - pa)
        s = (a + b + c) / 2.0
        area_sq = s * (s - a) * (s - b) * (s - c)
        if area_sq <= 0:
            continue
        if (a * b * c) / (4.0 * np.sqrt(area_sq)) < 1.0 / alpha:
            for i, j in [(0, 1), (1, 2), (2, 0)]:
                edges.add(tuple(sorted([simplex[i], simplex[j]])))
    lines = [((pts[i][0], pts[i][1]), (pts[j][0], pts[j][1])) for i, j in edges]
    polys = list(polygonize(MultiLineString(lines)))
    return unary_union(polys).buffer(buffer_dist, resolution=resolution) if polys else None


def smooth_hull_coords(hull, sigma=3, n_pts=1500):
    """Resample each boundary ring to n_pts equally spaced points, then Gaussian-smooth (wrap)."""
    if hull is None:
        return []
    polygons = [hull] if hull.geom_type == 'Polygon' else list(hull.geoms)
    out = []
    for p in polygons:
        x, y = np.array(p.exterior.coords.xy)
        seg = np.hypot(np.diff(x), np.diff(y)); t = np.concatenate([[0], np.cumsum(seg)])
        u = np.linspace(0, t[-1], n_pts, endpoint=False)
        xr, yr = np.interp(u, t, x), np.interp(u, t, y)
        out.append((gaussian_filter1d(xr, sigma, mode='wrap'), gaussian_filter1d(yr, sigma, mode='wrap')))
    return out


def compute_hull(xy):
    if len(xy) > SUBSAMPLE:
        xy = xy[np.random.default_rng(42).choice(len(xy), SUBSAMPLE, replace=False)]
    return smooth_hull_coords(concave_hull(xy, ALPHA, BUFFER, SMOOTH_RES), sigma=GAUSS_SIGMA)


def plot_hull(ax, hull, **kw):
    style = dict(linestyle=(0, (4, 3)), color='#444444', linewidth=0.9)
    style.update(kw)
    for x, y in hull:
        ax.plot(x, y, **style)


def hull_bounds(hull, pad=0.05):
    xs = np.concatenate([x for x, _ in hull]); ys = np.concatenate([y for _, y in hull])
    return xs.min() - pad, xs.max() + pad, ys.min() - pad, ys.max() + pad


# ---------------- data ----------------
def load(recompute=False):
    """DataFrame(sample, MDE_1, MDE_2, pluri_score) for all cells + hull coords."""
    if CACHE.exists() and not recompute:
        df = pd.read_parquet(CACHE)
        return df, compute_hull(df[['MDE_1', 'MDE_2']].to_numpy())
    import scanpy as sc
    adata = sc.read_h5ad(H5AD)
    missing = [g for g in PLURI_GENES if g not in adata.var_names]
    assert not missing, f'genes missing from var: {missing}'
    sc.tl.score_genes(adata, gene_list=PLURI_GENES, score_name=SCORE,
                      ctrl_size=50, n_bins=25, random_state=0)
    df = pd.DataFrame(adata.obsm['X_mde'], columns=['MDE_1', 'MDE_2'], index=adata.obs_names)
    df['sample'] = adata.obs['sample'].astype(str).to_numpy()
    df[SCORE] = adata.obs[SCORE].to_numpy()
    OUT.mkdir(exist_ok=True)
    df.to_parquet(CACHE)
    return df, compute_hull(adata.obsm['X_mde'])


def score_scatter(ax, sub, hull, bounds, vmin, vmax):
    """One panel: hull + cells coloured by score (high scores on top). Returns the mappable."""
    for spine in ax.spines.values():
        spine.set_visible(False)
    plot_hull(ax, hull)
    if POINT_ORDER == 'random':
        sub = sub.sample(frac=1, random_state=0)
    else:
        sub = sub.sort_values(SCORE)          # high scores on top
    m = ax.scatter(sub['MDE_1'], sub['MDE_2'], c=sub[SCORE], cmap=CMAP, vmin=vmin, vmax=vmax,
                   s=1.2, linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlim(bounds[0], bounds[1]); ax.set_ylim(bounds[2], bounds[3])
    ax.set_aspect('equal')
    return m


def curved_arrow(ax, p0, p1, rad, coords='data', lw=1.3, color='black', shrink=4, head=14, fig=None):
    """arc3 quadratic arrow. rad<0 bends towards the upper-left of the p0->p1 direction (up-then-right)."""
    fig = fig or ax.figure
    patch = FancyArrowPatch(p0, p1, connectionstyle=f'arc3,rad={rad}', arrowstyle='-|>',
                            mutation_scale=head, lw=lw, color=color, shrinkA=shrink, shrinkB=shrink,
                            transform=ax.transData if coords == 'data' else fig.transFigure,
                            clip_on=False, zorder=10)
    if coords == 'data':
        ax.add_patch(patch)
    else:
        fig.add_artist(patch)


def quarter_arrow(fig, start_frac, up=True, radius_in=0.28, lw=2.0, head_len_in=0.14, head_w_in=0.12, color='black'):
    """Quarter-circle shaft starting vertically at `start_frac` (figure fraction), ending horizontal, head pointing right.
    Drawn in inches (dpi_scale_trans) so the arc is a true circle and the head base is perpendicular to the shaft end."""
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch, Polygon
    W, H = fig.get_size_inches()
    sx, sy = start_frac[0] * W, start_frac[1] * H
    r = radius_in
    cx, cy = sx + r, sy                                   # circle centre to the right of the start point
    th = np.linspace(np.pi, np.pi / 2, 60) if up else np.linspace(np.pi, 3 * np.pi / 2, 60)
    xs, ys = cx + r * np.cos(th), cy + r * np.sin(th)     # from start (180°) to the top/bottom (90°/270°), tangent -> right
    shaft = PathPatch(Path(np.column_stack([xs, ys])), fill=False, lw=lw, color=color, capstyle='butt',
                      transform=fig.dpi_scale_trans, clip_on=False, zorder=10)
    tx, ty = xs[-1], ys[-1]
    head = Polygon([[tx, ty - head_w_in / 2], [tx, ty + head_w_in / 2], [tx + head_len_in, ty]], closed=True,
                   fc=color, ec=color, lw=0, transform=fig.dpi_scale_trans, clip_on=False, zorder=11)
    fig.add_artist(shaft); fig.add_artist(head)


# ---------------- Figure 1 (right part of fig-to-refine1.png) ----------------
# All positions are figure fractions measured from fig-to-refine1.png (right part, x >= 880 px of 1950).
PW, PH = 0.190, 0.277                        # panel width / height
COL_X = [0.220, 0.416, 0.617]                # three treatment columns
ROW_Y = {'12h': 0.510, '24h': 0.110}         # row bottoms
HPSC_BOX = [0.014, 0.320, 0.182, 0.277]      # hPSC panel


def figure1(df, hull, vmax, save_path):
    d = df[df['sample'].isin(FIG1_SAMPLES)]
    hp = d[d['sample'] == 'H1-hPSC.p6'][['MDE_1', 'MDE_2']].median().values
    bounds = hull_bounds(hull, pad=0.05)

    fig = plt.figure(figsize=(10.7, 6.1))
    fig.text(0.395, 0.930, SCORE_LABEL, fontsize=14, color=TITLE_COLOR, ha='right', va='center')
    fig.text(0.400, 0.930, '(single-cell RNA-seq)', fontsize=14, color='black', ha='left', va='center')

    layout = {'H1-0.5E': (COL_X[0], ROW_Y['12h']), 'H1-0.5ET': (COL_X[1], ROW_Y['12h']), 'H1-0.5E4T': (COL_X[2], ROW_Y['12h']),
              'H1-E': (COL_X[0], ROW_Y['24h']), 'H1-ET': (COL_X[1], ROW_Y['24h']), 'H1-E4T.1': (COL_X[2], ROW_Y['24h'])}
    axes = {'H1-hPSC.p6': fig.add_axes(HPSC_BOX)}
    for s_, (x, y) in layout.items():
        axes[s_] = fig.add_axes([x, y, PW, PH])
    m = None
    for s_, a in axes.items():
        m = score_scatter(a, d[d['sample'] == s_], hull, bounds, VMIN, vmax)
    # panel titles: one text block per panel, all centred on the same y (block centre), as in the reference
    for s_, (x, y) in layout.items():
        fig.text(x + PW / 2, y + PH + 0.055, PANEL_TITLE[s_], fontsize=13, ha='center', va='center',
                 linespacing=1.15)
    fig.text(HPSC_BOX[0] + 0.045, HPSC_BOX[1] + PH * 0.83, 'hPSC', fontsize=13, ha='left', va='center')

    # hPSC -> rows: exact quarter-circle shafts (up-then-right / down-then-right) + sharp triangular heads
    # measured on fig-to-refine1.png: starts 0.34 in above/below the hPSC panel centre, 0.07 in right of its edge
    x0 = HPSC_BOX[0] + HPSC_BOX[2] + 0.007
    yc = HPSC_BOX[1] + PH / 2
    yc -= 0.065                                   # pair shifted down so the upper head stays clear of the 12h panel
    quarter_arrow(fig, (x0, yc + 0.056), up=True)
    quarter_arrow(fig, (x0, yc - 0.056), up=False)

    # colourbar under the hPSC panel, labels on the left
    cax = fig.add_axes([HPSC_BOX[0] + 0.070, 0.215, 0.014, 0.100])   # tucked under the hull's left lobe
    cb = fig.colorbar(m, cax=cax)
    cb.outline.set_visible(False)
    cb.set_ticks([VMIN, vmax]); cb.set_ticklabels([f'{VMIN:g}', f'{vmax:g}'])
    cb.ax.tick_params(labelsize=12, length=0, pad=3)
    cb.ax.yaxis.set_ticks_position('left')

    # "Reversion" callout pointing at the hPSC position inside the 12h -> 4d panel
    ann = axes['H1-0.5E4T'].annotate('Reversion to\nhPSC-like state', xy=(hp[0] + 0.06, hp[1] + 0.02), xycoords='data',
                                     xytext=(COL_X[2] + PW + 0.012, 0.615), textcoords='figure fraction',
                                     fontsize=13.5, ha='left', va='center',
                                     bbox=dict(boxstyle='round,pad=0.45', fc='#d9d2e9', ec='none'),
                                     arrowprops=dict(arrowstyle='-|>', lw=1.3, color='black', shrinkA=2, shrinkB=2),
                                     annotation_clip=False)

    fig.canvas.draw()                                  # the callout's rounded box only exists after a draw
    extra = fig.get_default_bbox_extra_artists() + [ann, ann.get_bbox_patch()]
    fig.savefig(save_path, dpi=300, bbox_inches='tight', pad_inches=0.15, bbox_extra_artists=extra)
    plt.close(fig)
    print('Saved:', save_path)


def main():
    global POINT_ORDER, VMIN, VMAX_PERCENTILE
    if '--order' in sys.argv:
        POINT_ORDER = sys.argv[sys.argv.index('--order') + 1]
    if '--vmax-pct' in sys.argv:
        VMAX_PERCENTILE = float(sys.argv[sys.argv.index('--vmax-pct') + 1])
    if '--vmin' in sys.argv:
        VMIN = float(sys.argv[sys.argv.index('--vmin') + 1])
    tag = sys.argv[sys.argv.index('--tag') + 1] if '--tag' in sys.argv else ''
    df, hull = load(recompute='--recompute' in sys.argv)
    print(df.groupby('sample')[SCORE].agg(['mean', 'median', 'max']).round(3).to_string())
    vmax = float(np.round(np.percentile(df[SCORE], VMAX_PERCENTILE), 2))
    print(f'colorbar range: {VMIN} .. {vmax}')
    suffix = ('' if POINT_ORDER == 'score' else f'.{POINT_ORDER}order') + (f'.{tag}' if tag else '')
    figure1(df, hull, vmax, OUT / f'Figure1.PluriScore.MDEmap{suffix}.pdf')


if __name__ == '__main__':
    main()
