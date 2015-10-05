#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
import matplotlib
matplotlib.use('Agg')

import os.path as op
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from .plotting import plot_measures, plot_mosaic, plot_all, plot_fd

# matplotlib.rc('figure', figsize=(11.69, 8.27))  # for DINA4 size


def report_anatomical(in_csv, sc_split=False, split_files=True,
                      out_file='anatomical.pdf'):
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

    subject_list = pd.unique(df.subject.ravel())

    if not split_files:
        report = PdfPages(out_file)
    else:
        tpl, _ = op.splitext(op.basename(out_file))
        tpl = op.join(op.dirname(out_file), tpl) + '_%s.pdf'

    for ss in sessions:
        sesdf = df.loc[df['session'] == ss]
        scans = pd.unique(sesdf.scan.ravel())

        for subject in subject_list:
            if split_files:
                report = PdfPages(tpl % subject)

            subdf = sesdf.loc[sesdf['subject'] == subject]
            scans = pd.unique(subdf.scan.ravel())
            for sc in scans:
                data = subdf.loc[subdf['scan'] == sc].anat_data.tolist()
                fig = plot_mosaic(data[0], title=subject)
                report.savefig(fig, dpi=300)
                fig.clf()
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
                    fig = plot_all(sesdf, groups, subject=subject)
                    report.savefig(fig, dpi=300)
                    fig.clf()
                    fig = plot_measures(
                        sesdf, headers, subject=subject,
                        title='Report %s' % ss)
                    report.savefig(fig, dpi=300)
                    fig.clf()

            if split_files:
                report.close()
                print 'Written report for subject %s' % subject

    if not split_files:
        report.close()

    plt.close()
    return out_file


def report_func_temporal(in_csv, sc_split=False, split_files=True,
                         out_file='func_temporal.pdf'):
    df = pd.read_csv(in_csv)

    groups = [[['dvars'], ['gcor'], ['mean_tsnr'], ['mean_fd'],
               ['num_fd'], ['outlier'], ['perc_fd'], ['quality']]]
    headers = [v for gnames in groups for v in gnames]

    sessions = pd.unique(df.session.ravel())
    subject_list = pd.unique(df.subject.ravel())

    if not split_files:
        report = PdfPages(out_file)
    else:
        tpl, _ = op.splitext(op.basename(out_file))
        tpl = op.join(op.dirname(out_file), tpl) + '_%s.pdf'

    for ss in sessions:
        sesdf = df.loc[df['session'] == ss]
        scans = pd.unique(sesdf.scan.ravel())

        for subject in subject_list:
            if split_files:
                report = PdfPages(tpl % subject)

            subdf = sesdf.loc[sesdf['subject'] == subject]
            scans = pd.unique(subdf.scan.ravel())

            if 'fd_file' in sesdf.columns:
                for sc in scans:
                    data = subdf.loc[subdf['scan'] == sc].fd_file.tolist()[0]
                    mean_df = sesdf[['mean_df']].tolist()
                    fig = plot_fd(data, mean_df_dist=mean_df, title=subject)
                    report.savefig(fig, dpi=300)
                    fig.clf()

            if 'tsnr_file' in sesdf.columns:
                for sc in scans:
                    data = subdf.loc[subdf['scan'] == sc].tsnr_file.tolist()[0]
                    fig = plot_mosaic(data, title=subject)
                    report.savefig(fig, dpi=300)
                    fig.clf()

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
                    fig = plot_all(sesdf, groups, subject=subject)
                    report.savefig(fig, dpi=300)
                    fig.clf()
                    fig = plot_measures(
                        sesdf, headers, subject=subject,
                        title='Report %s' % ss)
                    report.savefig(fig, dpi=300)
                    fig.clf()

            if split_files:
                report.close()
                print 'Written report for subject %s' % subject

    if not split_files:
        report.close()

    plt.close()
    return out_file


def report_func_spatial(in_csv, sc_split=False, split_files=True,
                        out_file='func_spatial.pdf'):
    df = pd.read_csv(in_csv)

    groups = [[['dvars'], ['gcor'], ['mean_fd'],
               ['num_fd'], ['outlier'], ['perc_fd'], ['quality']]]
    headers = [v for gnames in groups for v in gnames]

    sessions = pd.unique(df.session.ravel())
    subject_list = pd.unique(df.subject.ravel())

    if not split_files:
        report = PdfPages(out_file)
    else:
        tpl, _ = op.splitext(op.basename(out_file))
        tpl = op.join(op.dirname(out_file), tpl) + '_%s.pdf'

    for ss in sessions:
        sesdf = df.loc[df['session'] == ss]
        scans = pd.unique(sesdf.scan.ravel())

        for subject in subject_list:
            if split_files:
                report = PdfPages(tpl % subject)

            subdf = sesdf.loc[sesdf['subject'] == subject]
            scans = pd.unique(subdf.scan.ravel())

            for sc in scans:
                data = subdf.loc[subdf['scan'] == sc].mean_epi.tolist()
                fig = plot_mosaic(data[0], title=subject)
                report.savefig(fig, dpi=300)
                fig.clf()

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
                    fig = plot_all(sesdf, groups, subject=subject)
                    report.savefig(fig, dpi=300)
                    fig.clf()
                    fig = plot_measures(
                        sesdf, headers, subject=subject,
                        title='Report %s' % ss)
                    report.savefig(fig, dpi=300)
                    fig.clf()

            if split_files:
                report.close()
                print 'Written report for subject %s' % subject

    if not split_files:
        report.close()

    plt.close()
    return out_file
