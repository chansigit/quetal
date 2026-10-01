"""Figure 2B (bottom): pluripotency gene-set score violins, Control vs JARID2-CRISPRi.

Replaces the NANOG violin panel of fig-to-refine2.png with the 8-gene pluripotency score
(same score as pluri_score_mde.py; read from figures/pluri_score.cache.parquet, run that script
first if the cache is missing). Colours and layout follow the reference panel.

    python 5d.pluripotency-score/pluri_score_violin_fig2.py
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans']
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent.parent      # workflow-250505/
THIS = Path(__file__).resolve().parent             # 5d.pluripotency-score/
OUT = HERE / 'figures'
CACHE = THIS / 'pluri_score.cache.parquet'
SCORE = 'pluri_score'

# (group header, control sample, JARID2 sample); colours from workflow-250505/4b + 5a
GROUPS = [
    ('hPSC',                            'Ctrl-hPSC.p1', 'JARID2-hPSC'),
    ('Ectoderm',                        'Ctrl-E.p6',    'JARID2-E'),
    ('24 h Ecto.\n→1d Plur.\nmedia',    'Ctrl-ET.p1',   'JARID2-ET'),
    ('24 h Ecto.\n→4d Plur.\nmedia',    'H1-E4T.3',     'JARID2-E4T'),
]
COLORS = {
    'Ctrl-hPSC.p1': '#f6c445', 'JARID2-hPSC': '#CE8E1B',
    'Ctrl-E.p6':    '#72a699', 'JARID2-E':    '#477343',
    'Ctrl-ET.p1':   '#3c69c3', 'JARID2-ET':   '#a1c8e7',
    'H1-E4T.3':     '#2e247e', 'JARID2-E4T':  '#78599d',
}
RED = '#e8231a'


def stars(p):
    return 'ns' if p >= 0.05 else '*' if p >= 0.01 else '**' if p >= 0.001 else '***'


def main():
    df = pd.read_parquet(CACHE)
    order = [s for _, c, j in GROUPS for s in (c, j)]
    d = df[df['sample'].isin(order)].copy()
    d['sample'] = pd.Categorical(d['sample'], order)

    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    sns.violinplot(data=d, x='sample', y=SCORE, order=order, hue='sample', hue_order=order, legend=False,
                   palette=COLORS, inner=None, cut=0, density_norm='width', width=0.8,
                   linewidth=0.8, linecolor='black', ax=ax)

    # median of each violin: white dot with black edge
    med = d.groupby('sample', observed=True)[SCORE].median()
    ax.scatter(range(len(order)), [med[s_] for s_ in order], s=22, c='white', edgecolors='black',
               linewidths=0.8, zorder=5)

    # axes style: open frame like the reference
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.spines['left'].set_linewidth(1.0); ax.spines['bottom'].set_linewidth(1.0)
    ax.set_xlabel('')
    ax.set_ylabel('Pluripotency score', fontsize=13)
    lo, hi = d[SCORE].min(), d[SCORE].max()
    ax.set_ylim(lo - 0.02, hi * 1.02)
    ax.set_yticks([0, round(hi, 1)]); ax.tick_params(axis='y', labelsize=12, length=3)
    ax.set_xlim(-0.6, len(order) - 0.4)

    # x tick labels: "Control" / "JARID2-CRISPRi" with JARID2 in red italic, rotated 45°
    ax.set_xticks(range(len(order))); ax.set_xticklabels([''] * len(order)); ax.tick_params(axis='x', length=0)
    for i, s in enumerate(order):
        if s.startswith('JARID2'):
            # two rotated lines sharing one anchor: 'JARID2' (red italic) and '-CRISPRi' one line-height below it,
            # i.e. offset perpendicular to the 45° baseline (down-right)
            common = dict(fontsize=12, rotation=45, ha='right', va='top', rotation_mode='anchor')
            ax.annotate('$\\it{JARID2}$', xy=(i, -0.03), xycoords=ax.get_xaxis_transform(), xytext=(0, 0),
                        textcoords='offset points', color=RED, **common)
            ax.annotate('-CRISPRi', xy=(i, -0.03), xycoords=ax.get_xaxis_transform(), xytext=(9.5, -9.5),
                        textcoords='offset points', color='black', **common)
        else:
            ax.text(i, -0.03, 'Control', fontsize=12, rotation=45, ha='right', va='top',
                    rotation_mode='anchor', transform=ax.get_xaxis_transform())

    # group headers with underline, significance (two-sided Mann-Whitney U) and delta median (JARID2 - Control)
    stats = []
    for k, (label, c, j) in enumerate(GROUPS):
        x0, x1 = 2 * k, 2 * k + 1
        a, b = d.loc[d['sample'] == c, SCORE], d.loc[d['sample'] == j, SCORE]
        p = mannwhitneyu(a, b, alternative='two-sided').pvalue   # two-sided; direction/effect size shown as delta median
        stats.append((label.replace('\n', ' '), c, j, len(a), len(b), a.median(), b.median(), b.median() - a.median(), p, stars(p)))
        xm = (x0 + x1) / 2
        ax.text(xm, 1.105, label, ha='center', va='bottom', fontsize=12, transform=ax.get_xaxis_transform(), linespacing=1.0)
        ax.plot([x0 - 0.35, x1 + 0.35], [1.095, 1.095], color='black', lw=1.0, transform=ax.get_xaxis_transform(), clip_on=False)
        ax.text(xm, 1.075, stars(p), ha='center', va='top', fontsize=12, transform=ax.get_xaxis_transform(),
                style='italic' if stars(p) == 'ns' else 'normal')
        ax.text(xm, 1.015, f'Δ = {b.median() - a.median():+.3f}', ha='center', va='top', fontsize=10,
                transform=ax.get_xaxis_transform())

    OUT.mkdir(exist_ok=True)
    out = OUT / 'Figure2.PluriScore.violin.pdf'
    fig.savefig(out, dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)
    print('Saved:', out)
    st = pd.DataFrame(stats, columns=['group', 'control', 'jarid2', 'n_control', 'n_jarid2', 'median_control', 'median_jarid2', 'delta_median', 'p_mannwhitney_twosided', 'label'])
    print(st.to_string(index=False))
    st.to_csv(THIS / 'Figure2.PluriScore.violin.stats.csv', index=False)


if __name__ == '__main__':
    main()
