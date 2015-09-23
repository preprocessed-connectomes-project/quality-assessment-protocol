#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
import seaborn as sns
import pandas as pd
from matplotlib import rc
import matplotlib.pyplot as plt

# rc('figure', figsize=(11.69, 8.27))  # for DINA4 size


def plot_all(df, groups, subject=None, out_file='testplot.pdf'):
    import matplotlib.gridspec as gridspec
    colnames = [v for gnames in groups for v in gnames]
    lengs = [len(el) for el in groups]
    ncols = np.sum(lengs)

    fig = plt.figure(figsize=[ncols, 5])
    gs = gridspec.GridSpec(1, len(groups), width_ratios=lengs)

    axes = []
    for i, snames in enumerate(groups):
        axes.append(plt.subplot(gs[i]))
        sns.violinplot(data=df[snames], ax=axes[-1])
        #df[snames].plot(kind='box', ax=axes[-1])

        if subject is not None:
            fig.suptitle(subject)
            out_file = subject + '.pdf'
            subdf = df.loc[df['subject'] == subject].copy()
            reps = len(subdf)
            for j, s in enumerate(snames):
                vals = np.atleast_1d(subdf[[s]]).reshape(-1).tolist()

                pos = np.array([j] * len(vals)) + 0.12 * \
                    (np.arange(0.0, len(vals)) - 0.5 * len(vals))
                axes[-1].plot(pos, vals, markersize=5, linestyle='None',
                              color='w', marker='o', markeredgecolor='k')

    fig.savefig(out_file, dpi=100, bbox_inches='tight')
    return fig


def plot_vline(cur_val, label, ax):
    ax.axvline(cur_val)
    ylim = ax.get_ylim()
    vloc = (ylim[0] + ylim[1]) / 2.0
    xlim = ax.get_xlim()
    pad = (xlim[0] + xlim[1]) / 100.0
    ax.text(cur_val - pad, vloc, label, color="blue", rotation=90,
            verticalalignment='center', horizontalalignment='right')


def plot_measures(df, measures, ncols=4, subject=None,
                  out_file='testplot.pdf'):
    import matplotlib.gridspec as gridspec
    nmeasures = len(measures)
    nrows = nmeasures // ncols
    if nmeasures % ncols > 0:
        nrows += 1

    fig = plt.figure(figsize=(3*ncols, 3*nrows))
    gs = gridspec.GridSpec(nrows, ncols)

    axes = []

    for i, mname in enumerate(measures):
        axes.append(plt.subplot(gs[i]))
        axes[-1].set_xlabel(mname)
        sns.distplot(df[[mname]], ax=axes[-1], color="b")

        if subject is not None:
            fig.suptitle(subject)

            subdf = df.loc[df['subject'] == subject].copy()
            sessions = np.atleast_1d(subdf[['session']]).reshape(-1).tolist()

            for ss in sessions:
                sesdf = subdf.loc[subdf['session'] == ss]
                scans = np.atleast_1d(sesdf[['scan']]).reshape(-1).tolist()

                for sc in scans:
                    scndf = subdf.loc[sesdf['scan'] == sc]
                    plot_vline(
                        scndf.iloc[0][mname], '%s_%s' % (ss, sc), axes[-1])

                out_file = subject + '_distplots.pdf'

    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    fig.savefig(out_file, dpi=100, bbox_inches='tight')
    return fig
