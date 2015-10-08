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

# matplotlib.rc('figure', figsize=(11.69, 8.27))  # for DINA4 size


def report_anatomical(in_csv, sc_split=False, split_files=True,
                      out_file='anatomical.pdf'):
    df = pd.read_csv(in_csv).fillna(0)
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
    df = pd.read_csv(in_csv).fillna(0)

    groups = [['dvars'], ['gcor'], ['m_tsnr'], ['mean_fd'],
              ['num_fd'], ['outlier'], ['perc_fd'], ['quality']]
    headers = [v for gnames in groups for v in gnames]

    sessions = sorted(pd.unique(df.session.ravel()))
    subject_list = pd.unique(df.subject.ravel())

    rootreport = PdfPages(out_file)

    for sub_id in subject_list:
        if split_files:
            tpl, _ = op.splitext(op.basename(out_file))
            tpl = op.join(op.dirname(out_file), tpl) + '_%s.pdf'
            report = PdfPages(tpl % sub_id)
        else:
            report = rootreport

        for ss in sessions:
            sesdf = df.copy().loc[df['session'] == ss]
            subdf = sesdf.copy().loc[sesdf['subject'] == sub_id]
            scans = pd.unique(subdf.scan.ravel())

            if 'tsnr_file' in sesdf.columns:
                for sc in scans:
                    subtitle = '(subject %s_%s_%s)' % (sub_id, ss, sc)
                    N = len(sesdf.m_tsnr)

                    data = subdf.loc[subdf['scan'] == sc].tsnr_file.tolist()[0]
                    fig = plot_mosaic(data, title='tSNR map ' + subtitle)
                    report.savefig(fig, dpi=300)
                    fig.clf()

                    mask = subdf.loc[subdf['scan'] == sc].mask_file.tolist()[0]
                    fig = plot_dist(
                        data, mask, 'tSNR inside the mask ' + subtitle,
                        sesdf.m_tsnr, ('Distribution of median tSNR of all'
                                       ' subjects (N=%d)') % N,
                        figsize=(8.3, 8.3))
                    report.savefig(fig, dpi=300)
                    fig.clf()

            if 'fd_file' in sesdf.columns:
                for sc in scans:
                    data = subdf.loc[subdf['scan'] == sc].fd_file.tolist()[0]
                    mean_fd = sesdf.mean_fd.tolist()
                    fig = plot_fd(data, mean_fd_dist=mean_fd)
                    report.savefig(fig, dpi=300)
                    fig.clf()

            if sc_split:
                for sc in scans:
                    subset = sesdf.loc[sesdf['scan'] == sc]
                    if len(subset.index) > 1:
                        fig = plot_measures(
                            subset, headers, subject=sub_id,
                            title='Report %s_%s' % (ss, sc),
                            figsize=(8.27, 4.06))
                        report.savefig(fig, dpi=300)
                        fig.clf()
            else:
                subtitle = '(subject %s_%s)' % (sub_id, ss)
                if len(sesdf.index) > 1:
                    fig = plot_all(sesdf, groups, subject=sub_id)
                    report.savefig(fig, dpi=300)
                    fig.clf()
                    fig = plot_measures(
                        sesdf, headers, subject=sub_id,
                        title='QC measures ' + subtitle,
                        figsize=(8.27, 4.06))
                    report.savefig(fig, dpi=300)
                    fig.clf()

        if split_files:
            report.close()
            fname = tpl % sub_id
            print 'Written report for subject (%s)' % fname

    # Group plots
    for ss in sessions:
        sesdf = df.copy().loc[df['session'] == ss]

        if len(sesdf.index) > 1:
            fig = plot_measures(
                sesdf, headers, title='QC measures ' + ss,
                figsize=(8.27, 4.06))
            rootreport.savefig(fig, dpi=300)
            fig.clf()
            fig = plot_all(sesdf, groups)
            rootreport.savefig(fig, dpi=300)
            fig.clf()

    rootreport.close()
    plt.close()
    return out_file


def report_func_spatial(in_csv, sc_split=False, split_files=True,
                        out_file='func_spatial.pdf'):
    df = pd.read_csv(in_csv).fillna(0)

    groups = [['bg_size', 'fg_size'],
              ['bg_mean', 'fg_mean'],
              ['bg_std', 'fg_std'],
              ['efc'],
              ['fber'],
              ['fwhm', 'fwhm_x', 'fwhm_y', 'fwhm_z'],
              ['ghost_x'],
              ['snr']]

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
