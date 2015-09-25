#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import math
import nibabel as nb
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.backends.backend_pdf import FigureCanvasPdf as FigureCanvas

# matplotlib.rc('figure', figsize=(11.69, 8.27))  # for DINA4 size


def _calc_rows_columns(ratio, n_images):
    rows = 1
    for _ in range(100):
        columns = math.floor(ratio * rows)
        total = rows * columns
        if total > n_images:
            break

        columns = math.ceil(ratio * rows)
        total = rows * columns
        if total > n_images:
            break
        rows += 1
    return rows, columns


def plot_vline(cur_val, label, ax):
    ax.axvline(cur_val)
    ylim = ax.get_ylim()
    vloc = (ylim[0] + ylim[1]) / 2.0
    xlim = ax.get_xlim()
    pad = (xlim[0] + xlim[1]) / 100.0
    ax.text(cur_val - pad, vloc, label, color="blue", rotation=90,
            verticalalignment='center', horizontalalignment='right')


def plot_measures(df, measures, ncols=4, title='Group level report',
                  subject=None, out_file='testplot.pdf'):
    import matplotlib.gridspec as gridspec
    nmeasures = len(measures)
    nrows = nmeasures // ncols
    if nmeasures % ncols > 0:
        nrows += 1

    fig = plt.figure(figsize=(3 * ncols, 3 * nrows))
    gs = gridspec.GridSpec(nrows, ncols)

    axes = []

    for i, mname in enumerate(measures):
        axes.append(plt.subplot(gs[i]))
        axes[-1].set_xlabel(mname)
        sns.distplot(df[[mname]], ax=axes[-1], color="b")

        if subject is not None:
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

    if subject is not None:
        fig.suptitle(subject)
    else:
        fig.suptitle(title)

    # plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    return fig


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
        # df[snames].plot(kind='box', ax=axes[-1])

        if subject is not None:
            fig.suptitle(subject)
            out_file = subject + '.pdf'
            subdf = df.loc[df['subject'] == subject].copy()
            reps = len(subdf)
            for j, s in enumerate(snames):
                vals = np.atleast_1d(subdf[[s]]).reshape(-1).tolist()

                pos = [j]
                if len(vals) > 1:
                    pos = np.array([j] * len(vals)) + 0.12 * \
                        (np.arange(0.0, len(vals)) - 0.5 * len(vals))

                axes[-1].plot(pos, vals, markersize=5, linestyle='None',
                              color='w', marker='o', markeredgecolor='k')
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    return fig


def plot_mosaic(nifti_file, title=None, overlay_mask=None,
                figsize=(11.7, 8.3)):
    if isinstance(nifti_file, str):
        nii = nb.load(nifti_file)
        mean_data = nii.get_data()
    else:
        mean_data = nifti_file

    n_images = mean_data.shape[2]
    row, col = _calc_rows_columns(figsize[0] / figsize[1], n_images)

    if overlay_mask:
        overlay_data = nb.load(overlay_mask).get_data()

    # create figures
    fig = Figure(figsize=figsize)
    FigureCanvas(fig)

    fig.subplots_adjust(top=0.85)
    for image in (range(n_images)):
        ax = fig.add_subplot(row, col, image + 1)
        data_mask = np.logical_not(np.isnan(mean_data))
        if overlay_mask:
            ax.set_rasterized(True)

        ax.imshow(np.fliplr(mean_data[:, :, image].T), vmin=np.percentile(
            mean_data[data_mask], 0.5),
            vmax=np.percentile(mean_data[data_mask], 99.5),
            cmap=cm.Greys_r, interpolation='nearest', origin='lower')
        if overlay_mask:
            cmap = cm.Reds  # @UndefinedVariable
            cmap._init()
            alphas = np.linspace(0, 0.75, cmap.N + 3)
            cmap._lut[:, -1] = alphas
            ax.imshow(np.fliplr(overlay_data[:, :, image].T), vmin=0, vmax=1,
                      cmap=cmap, interpolation='nearest', origin='lower')

        ax.axis('off')

    fig.subplots_adjust(
        left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=0.01, hspace=0.1)

    if not title:
        _, title = os.path.split(nifti_file)
        title += " (last modified: %s)" % time.ctime(
            os.path.getmtime(nifti_file))
    fig.suptitle(title, fontsize='10')

    return fig


def report_anatomical(in_csv, subject=None, sc_split=False,
                      out_file='anatomical.pdf'):
    import numpy as np
    import pandas as pd
    import math
    import nibabel as nb
    import seaborn as sns
    from matplotlib import rc
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.backends.backend_pdf import FigureCanvasPdf as FigureCanvas
    report = PdfPages(out_file)
    df = pd.read_csv(in_csv)
    sessions = pd.unique(df.session.ravel())

    groups = [['bg_size', 'fg_size'],
              ['bg_mean', 'fg_mean'],
              ['bg_std', 'fg_std'],
              ['csf_size', 'gm_size', 'wm_size'],
              ['csf_mean', 'gm_mean', 'wm_mean'],
              ['csf_std', 'gm_std', 'wm_std'],
              ['cnr'],
              ['efc'],
              ['fber'],
              ['fwhm', 'fwhm_x', 'fwhm_y', 'fwhm_z'],
              ['qi1'],
              ['snr']]
    headers = [v for gnames in groups for v in gnames]

    for ss in sessions:
        sesdf = df.loc[df['session'] == ss]
        scans = pd.unique(sesdf.scan.ravel())

        if subject is not None:
            for sc in scans:
                data = sesdf.loc[sesdf['scan'] == sc]['anat_data']
                plot_mosaic(data, title=subject)

        if sc_split:
            for sc in scans:
                subset = sesdf.loc[sesdf['scan'] == sc]

                if len(subset.index) > 1:
                    fig = plot_measures(
                        subset, headers, subject=subject,
                        title='Report %s_%s' % (ss, sc))
                    report.savefig(fig, dpi=300)
                    fig.clf()
        else:
            if len(sesdf.index) > 1:
                fig = plot_measures(
                    sesdf, headers, subject=subject,
                    title='Report %s' % ss)
                report.savefig(fig, dpi=300)
                fig.clf()

                fig = plot_all(sesdf, groups, subject=subject)
                report.savefig(fig, dpi=300)
                fig.clf()

    report.close()
    plt.close()
    return out_file


def report_functional(in_csv, subject=None, sc_split=False,
                      out_file='functional.pdf'):
    import numpy as np
    import pandas as pd
    import math
    import nibabel as nb
    import seaborn as sns
    from matplotlib import rc
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.backends.backend_pdf import FigureCanvasPdf as FigureCanvas
    report = PdfPages(out_file)
    df = pd.read_csv(in_csv)
    sessions = pd.unique(df.session.ravel())

    groups = [[['dvars'], ['gcor'], ['mean_fd'],
               ['num_fd'], ['outlier'], ['perc_fd'], ['quality']]]
    headers = [v for gnames in groups for v in gnames]

    for ss in sessions:
        sesdf = df.loc[df['session'] == ss]

        if sc_split:
            scans = pd.unique(sesdf.scan.ravel())

            for sc in scans:
                subset = sesdf.loc[sesdf['scan'] == sc]

                if len(subset.index) > 1:
                    fig = plot_measures(
                        subset, headers, subject=subject,
                        title='Report %s_%s' % (ss, sc))
                    report.savefig(fig, dpi=300)
                    fig.clf()
        else:
            if len(sesdf.index) > 1:
                fig = plot_measures(
                    sesdf, headers, subject=subject,
                    title='Report %s' % ss)
                report.savefig(fig, dpi=300)
                fig.clf()

                fig = plot_all(sesdf, groups, subject=subject)
                report.savefig(fig, dpi=300)
                fig.clf()

    report.close()
    plt.close()
    return out_file
