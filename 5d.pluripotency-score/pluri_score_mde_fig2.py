"""Figure 2B (top): Control vs JARID2-CRISPRi MDE panels coloured by the pluripotency score.

Same two-panel layout, trajectory arrows, cluster labels and MDE axis indicator as
5a (Figure2.Samples.MDEmap.pdf); only the point colour changes from cell
type to the 8-gene pluripotency score (5d.pluripotency-score/pluri_score.cache.parquet, written by
pluri_score_mde.py). No hull. Colour range: 0.05 .. 99.5th percentile of these 8 samples.

    python 5d.pluripotency-score/pluri_score_mde_fig2.py
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans']
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea

HERE = Path(__file__).resolve().parent.parent      # workflow-250505/
THIS = Path(__file__).resolve().parent             # 5d.pluripotency-score/
OUT = HERE / 'figures'
CACHE = THIS / 'pluri_score.cache.parquet'
SCORE = 'pluri_score'

CTRL = ['Ctrl-hPSC.p1', 'Ctrl-E.p6', 'Ctrl-ET.p1', 'H1-E4T.3']
JARID2 = ['JARID2-hPSC', 'JARID2-E', 'JARID2-ET', 'JARID2-E4T']

CMAP = mcolors.LinearSegmentedColormap.from_list('custom_cmap', ['#d1cfd4', '#9577e5', '#1206f5'])
VMIN = 0.05
VMAX_PERCENTILE = 99.5
RED = '#e8231a'

def anchors(sub, hpsc, ecto, et):
    """Data-space anchor points: medians of the hPSC and ectoderm samples, and the tip of the left arm."""
    hp = sub.loc[sub['sample'] == hpsc, ['MDE_1', 'MDE_2']].median().to_numpy()
    ec = sub.loc[sub['sample'] == ecto, ['MDE_1', 'MDE_2']].median().to_numpy()
    arm = sub.loc[sub['sample'] == et]
    other = arm.nsmallest(max(20, len(arm) // 20), 'MDE_1')[['MDE_1', 'MDE_2']].mean().to_numpy()
    return hp, ec, other


def arrow(ax, p0, p1, rad):
    ax.annotate('', xy=tuple(p1), xytext=tuple(p0), xycoords='data', textcoords='data',
                arrowprops=dict(arrowstyle='->', mutation_scale=25, lw=1.6, color='black',
                                connectionstyle=f'arc3,rad={rad}', shrinkA=0, shrinkB=0))


def panel(ax, sub, vmax, hpsc, ecto, et, reversion):
    sub = sub.sort_values(SCORE)                       # high scores on top
    m = ax.scatter(sub['MDE_1'], sub['MDE_2'], c=sub[SCORE], cmap=CMAP, vmin=VMIN, vmax=vmax,
                   s=5, linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlabel(''); ax.set_ylabel('')

    hp, ec, other = anchors(sub, hpsc, ecto, et)
    lab = dict(fontsize=20, va='center')
    ax.text(hp[0] + 0.22, hp[1] - 0.02, 'hPSC', ha='left', **lab)
    ax.text(ec[0] + 0.20, ec[1] + 0.02, 'Ectoderm', ha='left', **lab)
    ax.text(other[0] - 0.02, other[1] + 0.22, 'Other', ha='left', **lab)
    # hPSC -> Ectoderm: up along the right, bowing right (as in the reference)
    arrow(ax, hp + [0.08, 0.22], ec + [0.06, -0.25], rad=-0.2)
    if reversion:   # JARID2: Ectoderm -> hPSC, a wide arc on the left of the up-arrow, head pointing down into hPSC
        arrow(ax, ec + [-0.18, -0.20], hp + [-0.16, 0.20], rad=0.45)
    else:           # Control: Ectoderm -> Other, arc bowing down-left, head pointing left
        arrow(ax, ec + [-0.15, -0.15], other + [0.25, 0.05], rad=-0.3)   # bows below the chord, as in the reference
    return m


def main():
    df = pd.read_parquet(CACHE)
    d = df[df['sample'].isin(CTRL + JARID2)]
    vmax = float(np.round(np.percentile(d[SCORE], VMAX_PERCENTILE), 2))
    print(f'colorbar range: {VMIN} .. {vmax}  (99.5th pct of the 8 Figure 2 samples)')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7), sharex=True, sharey=True)
    m = panel(ax1, d[d['sample'].isin(CTRL)], vmax, 'Ctrl-hPSC.p1', 'Ctrl-E.p6', 'Ctrl-ET.p1', reversion=False)
    panel(ax2, d[d['sample'].isin(JARID2)], vmax, 'JARID2-hPSC', 'JARID2-E', 'JARID2-ET', reversion=True)

    ax1.set_title('Control', fontsize=22)
    t1 = TextArea('JARID2', textprops=dict(color=RED, fontstyle='italic', fontsize=22))
    t2 = TextArea('-CRISPRi', textprops=dict(color='black', fontsize=22))
    ax2.add_artist(AnchoredOffsetbox(loc='lower center', child=HPacker(children=[t1, t2], pad=0, sep=0),
                                     bbox_to_anchor=(0.5, 1.0), bbox_transform=ax2.transAxes, frameon=False))

    # MDE axis indicator (Control panel, lower left)
    origin = (0.03, 0.05)
    for tip in [(0.14, 0.05), (0.03, 0.19)]:   # both arrows start exactly at the origin (no shrink), so the tails join
        ax1.annotate('', xy=tip, xytext=origin, xycoords='axes fraction',
                     arrowprops=dict(arrowstyle='-|>', color='black', lw=1.2, shrinkA=0, shrinkB=0,
                                     mutation_scale=14, joinstyle='miter', capstyle='projecting'))
    ax1.text(0.085, 0.01, 'MDE 1', transform=ax1.transAxes, fontsize=14, ha='center', va='top')
    ax1.text(-0.01, 0.12, 'MDE 2', transform=ax1.transAxes, fontsize=14, ha='right', va='center', rotation=90)

    # shared colourbar on the right
    fig.subplots_adjust(left=0.03, right=0.90, wspace=0.08)
    cb = fig.colorbar(m, cax=fig.add_axes([0.915, 0.35, 0.012, 0.30]))
    cb.outline.set_visible(False)
    cb.set_ticks([VMIN, vmax]); cb.set_ticklabels([f'{VMIN:g}', f'{vmax:g}'])
    cb.ax.tick_params(labelsize=16, length=0)
    cb.set_label('Pluripotency score', fontsize=16)

    OUT.mkdir(exist_ok=True)
    out = OUT / 'Figure2.PluriScore.MDEmap.pdf'
    fig.savefig(out, dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)
    print('Saved:', out)


if __name__ == '__main__':
    main()
