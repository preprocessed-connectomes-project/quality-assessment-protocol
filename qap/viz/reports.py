#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import os.path as op
import pandas as pd
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from .plotting import (plot_measures, plot_mosaic, plot_all,
                       plot_fd, plot_dist)

# matplotlib.rc('figure', figsize=(11.69, 8.27))  # for DINA4 size


def get_documentation(doc_type, out_file):
    import codecs
    import StringIO
    from xhtml2pdf import pisa
    # open output file for writing (truncated binary)
    result = open(out_file, "w+b")

    html_dir = op.abspath(
        op.join(op.dirname(__file__), 'html', '%s.html' % doc_type))

    with codecs.open(html_dir, mode='r', encoding='utf-8') as f:
        html = f.read()

    # convert HTML to PDF
    status = pisa.pisaDocument(html, result, encoding='UTF-8')
    result.close()

    # return True on success and False on errors
    return status.err


def _read_csv(in_csv):
    df = pd.read_csv(in_csv, dtype={'subject': str}).fillna(
        0).drop_duplicates().sort(['subject', 'session'])
    return df


def concat_pdf(in_files, out_file='concatenated.pdf'):
    """
    Concatenate PDF list (http://stackoverflow.com/a/3444735)
    """
    from PyPDF2 import PdfFileWriter, PdfFileReader
    outpdf = PdfFileWriter()

    for in_file in in_files:
        inpdf = PdfFileReader(file(in_file, 'rb'))
        for p in range(inpdf.numPages):
            outpdf.addPage(inpdf.getPage(p))
    outpdf.write(file(out_file, 'wb'))
    return out_file


def _write_report(df, groups, sub_id=None, sc_split=False, condensed=True,
                  out_file='report.pdf'):
    columns = df.columns.ravel()
    headers = []
    for g in groups:
        rem = []
        for h in g:
            if h not in columns:
                rem.append(h)
            else:
                headers.append(h)
        for r in rem:
            g.remove(r)

    report = PdfPages(out_file)
    sessions = sorted(pd.unique(df.session.ravel()))
    for ss in sessions:
        sesdf = df.copy().loc[df['session'] == ss]
        scans = pd.unique(sesdf.scan.ravel())
        if sc_split:
            for sc in scans:
                subset = sesdf.loc[sesdf['scan'] == sc]
                if len(subset.index) > 1:
                    if sub_id is None:
                        subtitle = '(%s_%s)' % (ss, sc)
                    else:
                        subtitle = '(subject %s_%s_%s)' % (sub_id, ss, sc)
                    if condensed:
                        fig = plot_all(sesdf, groups, subject=sub_id,
                                       title='QC measures ' + subtitle)
                    else:
                        fig = plot_measures(
                            sesdf, headers, subject=sub_id,
                            title='QC measures ' + subtitle)
                    report.savefig(fig, dpi=300)
                    fig.clf()
        else:
            if len(sesdf.index) > 1:
                if sub_id is None:
                    subtitle = '(%s)' % (ss)
                else:
                    subtitle = '(subject %s_%s)' % (sub_id, ss)
                if condensed:
                    fig = plot_all(sesdf, groups, subject=sub_id,
                                   title='QC measures ' + subtitle)
                else:
                    fig = plot_measures(
                        sesdf, headers, subject=sub_id,
                        title='QC measures ' + subtitle)
                report.savefig(fig, dpi=300)
                fig.clf()

    report.close()
    plt.close()
    # print 'Written report file %s' % out_file
    return out_file


def _write_all_reports(df, groups, sc_split=False, condensed=True,
                       out_file='report.pdf'):

    outlist = []
    _write_report(
        df, groups, sc_split=sc_split, condensed=condensed, out_file=out_file)

    subject_list = sorted(pd.unique(df.subject.ravel()))
    for sub_id in subject_list:
        tpl, _ = op.splitext(op.basename(out_file))
        tpl = op.join(op.dirname(out_file), tpl) + '_%s.pdf'
        outlist.append(_write_report(
            df, groups, sub_id=sub_id, sc_split=sc_split, condensed=condensed,
            out_file=tpl % sub_id))
    return out_file, outlist


def all_anatomical(in_csv, sc_split=False, condensed=True,
                   out_file='anatomical.pdf'):
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
        _read_csv(in_csv), groups, sc_split=sc_split,
        condensed=condensed, out_file=out_file)


def all_func_temporal(in_csv, sc_split=False, condensed=True,
                      out_file='func_temporal.pdf'):
    groups = [['dvars'], ['gcor'], ['m_tsnr'], ['mean_fd'],
              ['num_fd'], ['outlier'], ['perc_fd'], ['quality']]
    return _write_all_reports(
        _read_csv(in_csv), groups, sc_split=sc_split,
        condensed=condensed, out_file=out_file)


def all_func_spatial(in_csv, sc_split=False, condensed=False,
                     out_file='func_spatial.pdf'):
    groups = [['bg_size', 'fg_size'],
              ['bg_mean', 'fg_mean'],
              ['bg_std', 'fg_std'],
              ['efc'],
              ['fber'],
              ['fwhm', 'fwhm_x', 'fwhm_y', 'fwhm_z'],
              ['ghost_%s' % a for a in ['x', 'y', 'z']],
              ['snr']]
    return _write_all_reports(
        _read_csv(in_csv), groups, sc_split=sc_split,
        condensed=condensed, out_file=out_file)


def report_anatomical(in_csv, subject=None, sc_split=False, condensed=True,
                      out_file='anatomical.pdf'):
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
    return _write_report(
        _read_csv(in_csv), groups, sub_id=subject,
        sc_split=sc_split, condensed=condensed, out_file=out_file)


def report_func_temporal(in_csv, subject=None, sc_split=False, condensed=True,
                         out_file='func_temporal.pdf'):
    groups = [['dvars'], ['gcor'], ['m_tsnr'], ['mean_fd'],
              ['num_fd'], ['outlier'], ['perc_fd'], ['quality']]
    return _write_report(
        _read_csv(in_csv), groups, sub_id=subject,
        sc_split=sc_split, condensed=condensed, out_file=out_file)


def report_func_spatial(in_csv, subject=None, sc_split=False, condensed=True,
                        out_file='func_spatial.pdf'):
    groups = [['bg_size', 'fg_size'],
              ['bg_mean', 'fg_mean'],
              ['bg_std', 'fg_std'],
              ['efc'],
              ['fber'],
              ['fwhm', 'fwhm_x', 'fwhm_y', 'fwhm_z'],
              ['ghost_%s' % a for a in ['x', 'y', 'z']],
              ['snr']]
    return _write_report(
        _read_csv(in_csv), groups, sub_id=subject,
        sc_split=sc_split, condensed=condensed, out_file=out_file)

