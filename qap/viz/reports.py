#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
import matplotlib
matplotlib.use('Agg')

import os.path as op
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from .plotting import (plot_measures, plot_mosaic, plot_all,
                       plot_fd, plot_dist)
from pyPdf import PdfFileWriter, PdfFileReader

# matplotlib.rc('figure', figsize=(11.69, 8.27))  # for DINA4 size


def concat_pdf(in_files, out_file='concatenated.pdf'):
    outpdf = PdfFileWriter()

    for in_file in in_files:
        inpdf = PdfFileReader(file(in_file, 'rb'))
        for p in range(inpdf.numPages):
            outpdf.addPage(inpdf.getPage(p))
    outpdf.write(file(out_file, 'wb'))
    return out_file


def _write_report(df, groups, sub_id=None, sc_split=False, condensed=True,
                  out_file='report.pdf'):
    headers = [v for gnames in groups for v in gnames]
    sessions = sorted(pd.unique(df.session.ravel()))

    report = PdfPages(out_file)
    for ss in sessions:
        sesdf = df.copy().loc[df['session'] == ss]
        subdf = sesdf.copy().loc[sesdf['subject'] == sub_id]
        scans = pd.unique(subdf.scan.ravel())
        if sc_split:
            for sc in scans:
                subset = sesdf.loc[sesdf['scan'] == sc]
                if len(subset.index) > 1:
                    subtitle = '(subject %s_%s_%s)' % (sub_id, ss, sc)
                    if condensed:
                        fig = plot_measures(
                            sesdf, headers, subject=sub_id,
                            title='QC measures ' + subtitle)
                    else:
                        fig = plot_all(sesdf, groups, subject=sub_id,
                                       title='QC measures ' + subtitle)
                    report.savefig(fig, dpi=300)
                    fig.clf()
        else:
            if len(sesdf.index) > 1:
                subtitle = '(subject %s_%s)' % (sub_id, ss)
                if condensed:
                    fig = plot_measures(
                        sesdf, headers, subject=sub_id,
                        title='QC measures ' + subtitle)
                else:
                    fig = plot_all(sesdf, groups, subject=sub_id,
                                   title='QC measures ' + subtitle)
                report.savefig(fig, dpi=300)
                fig.clf()

    report.close()
    plt.close()
    print 'Written report file %s' % out_file
    return out_file


def _write_all_reports(df, groups, sc_split=False, condensed=True,
                       out_file='report.pdf'):

    outlist = []
    _write_report(
        df, groups, sub_id=sub_id, sc_split=sc_split, condensed=condensed,
        out_file=out_file)

    for sub_id in subject_list:
        tpl, _ = op.splitext(op.basename(out_file))
        tpl = op.join(op.dirname(out_file), tpl) + '_%s.pdf'
        outlist.append(_write_report(
            df, groups, sub_id=sub_id, sc_split=sc_split, condensed=condensed,
            out_file=tpl))
    return out_file, outlist


def report_anatomical(in_csv, sc_split=False, condensed=True,
                      split_files=True, out_file='anatomical.pdf'):
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
    return _write_all_reports(
        pd.read_csv(in_csv).fillna(0), groups, sc_split=sc_split,
        condensed=condensed, split_files=split_files, out_file=out_file)


def report_func_temporal(in_csv, sc_split=False, split_files=True,
                         condensed=True, out_file='func_temporal.pdf'):
    groups = [['dvars'], ['gcor'], ['m_tsnr'], ['mean_fd'],
              ['num_fd'], ['outlier'], ['perc_fd'], ['quality']]
    return _write_all_reports(
        pd.read_csv(in_csv).fillna(0), groups, sc_split=sc_split,
        condensed=condensed, split_files=split_files, out_file=out_file)


def report_func_spatial(in_csv, sc_split=False, split_files=True,
                        condensed=False, out_file='func_spatial.pdf'):
    groups = [['bg_size', 'fg_size'],
              ['bg_mean', 'fg_mean'],
              ['bg_std', 'fg_std'],
              ['efc'],
              ['fber'],
              ['fwhm', 'fwhm_x', 'fwhm_y', 'fwhm_z'],
              ['ghost_x'],
              ['snr']]
    return _write_all_reports(
        pd.read_csv(in_csv).fillna(0), groups, sc_split=sc_split,
        condensed=condensed, split_files=split_files, out_file=out_file)
