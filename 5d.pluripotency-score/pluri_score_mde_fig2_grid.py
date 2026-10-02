"""Figure 2 (supplement): per-sample pluripotency-score MDE panels, Control (top row) vs
JARID2-CRISPRi (bottom row), one column per condition, in the style of the Figure 1G panels
(dashed concave-hull outline of all 27 samples, same colormap). Colour range 0.05 .. 99.5th
percentile of the 8 samples shown (as in Figure2.PluriScore.MDEmap.pdf).

    python 5d.pluripotency-score/pluri_score_mde_fig2_grid.py
"""
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pluri_score_mde as base   # hull, colormap, score_scatter, cache loader

HERE = Path(__file__).resolve().parent.parent      # workflow-250505/
OUT = HERE / 'figures'
RED = '#e8231a'

COLS = ['hPSC', '24h Ectoderm', '24h Ecto.→\n1d Plur. media', '24h Ecto.→\n4d Plur. media']
ROWS = [('Control', ['Ctrl-hPSC.p1', 'Ctrl-E.p6', 'Ctrl-ET.p1', 'H1-E4T.3']),
        ('JARID2', ['JARID2-hPSC', 'JARID2-E', 'JARID2-ET', 'JARID2-E4T'])]


def main():
    df, hull = base.load()
    samples = [s for _, r in ROWS for s in r]
    d = df[df['sample'].isin(samples)]
    vmin = 0.05
    # upper bound from the two 4d samples (95th pct), so that high-scoring reverted cells saturate
    d4 = d[d['sample'].isin(['H1-E4T.3', 'JARID2-E4T'])]
    VMAX_PCT = float(sys.argv[sys.argv.index('--vmax-pct') + 1]) if '--vmax-pct' in sys.argv else 99
    vmax = float(np.round(np.percentile(d4[base.SCORE], VMAX_PCT), 2))
    tag = sys.argv[sys.argv.index('--tag') + 1] if '--tag' in sys.argv else ''
    THR = 0.1   # per-panel annotation: fraction of cells above this score
    print(f'colorbar range: {vmin} .. {vmax}  ({VMAX_PCT:g}th pct of the two 4d samples)')
    bounds = base.hull_bounds(hull, pad=0.05)

    fig, axes = plt.subplots(2, 4, figsize=(11.2, 6.0))
    m = None
    for (row_label, row_samples), axrow in zip(ROWS, axes):
        for s_, ax in zip(row_samples, axrow):
            sub = d[d['sample'] == s_]
            m = base.score_scatter(ax, sub, hull, bounds, vmin, vmax)
            frac = (sub[base.SCORE] > THR).mean()
            ax.text(0.02, 0.97, f'score > {THR:g}:\n{frac:.0%} of cells', transform=ax.transAxes,
                    fontsize=10, ha='left', va='top', linespacing=1.1)
    for ax, title in zip(axes[0], COLS):
        ax.set_title(title, fontsize=13, pad=6, linespacing=1.15)
    # row labels at the left; JARID2 in red italic + "-CRISPRi"
    axes[0, 0].text(-0.08, 0.5, 'Control', transform=axes[0, 0].transAxes, fontsize=14,
                    ha='right', va='center', rotation=90)
    t1 = TextArea('JARID2', textprops=dict(color=RED, fontstyle='italic', fontsize=14, rotation=90))
    t2 = TextArea('-CRISPRi', textprops=dict(color='black', fontsize=14, rotation=90))
    from matplotlib.offsetbox import VPacker
    box = VPacker(children=[t2, t1], pad=0, sep=0, align='center')   # bottom-to-top reading order when rotated
    axes[1, 0].add_artist(AnchoredOffsetbox(loc='center right', child=box, frameon=False,
                                            bbox_to_anchor=(-0.02, 0.5), bbox_transform=axes[1, 0].transAxes))

    fig.subplots_adjust(left=0.07, right=0.90, top=0.90, bottom=0.03, wspace=0.08, hspace=0.18)
    cax = fig.add_axes([0.935, 0.42, 0.010, 0.16])
    cb = fig.colorbar(m, cax=cax)
    cb.outline.set_visible(False)
    cb.set_ticks([vmin, vmax]); cb.set_ticklabels([f'{vmin:g}', f'{vmax:g}'])
    cb.ax.tick_params(labelsize=11, length=0)
    cb.ax.yaxis.set_ticks_position('right')
    cb.set_label('Pluripotency score', fontsize=12)
    cb.ax.yaxis.set_label_position('left')

    out = OUT / f'Figure2.PluriScore.MDEmap.grid{"." + tag if tag else ""}.pdf'
    fig.savefig(out, dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)
    print('Saved:', out)


if __name__ == '__main__':
    main()
