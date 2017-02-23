#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# @Author: oesteban
# @Date:   2015-11-13 07:54:38

import sys
import os
import os.path as op
import pandas as pd
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from .plotting import (plot_measures, plot_mosaic, plot_all,
                       plot_fd, plot_dist)


def workflow_report(in_csv, qap_type, run_name, out_dir=None, out_file=None,
                    full_reports=False):
    """Generate a PDF report of a QAP run.

    :type in_csv: str
    :param in_csv: The filepath of the QAP output CSV file to create a report
                   for.
    :type qap_type: str
    :param qap_type: The type of QAP set of measures ("anatomical spatial",
                     etc.).
    :type run_name: str
    :param run_name: The name of the pipeline run.
    :param res_dict:
    :type out_dir: str
    :param out_dir: The output directory for the reports.
    :type out_file: str
    :param out_file: The filename of the PDF report.
    :type full_reports: bool
    :param full_reports: Whether or not to produce the individual-level
                         reports as well.
    :rtype: dict
    :return: A dictionary with information about the report generation.
    """

    if out_dir is None:
        out_dir = os.getcwd()

    if out_file is None:
        out_file = op.join(
            out_dir, qap_type + '_%s.pdf')

    # Read csv file, sort and drop duplicates
    df = pd.read_csv(in_csv, dtype={'Participant': str}).sort(
        columns=['Participant', 'Session', 'Series'])

    try:
        df.drop_duplicates(['Participant', 'Session', 'Series'], keep='last',
                           inplace=True)
    except TypeError:
        df.drop_duplicates(['Participant', 'Session', 'Series'],
                           take_last=True, inplace=True)

    df["Participant"] = df["Participant"].astype(str)
    df["Session"] = df["Session"].astype(str)
    df["Series"] = df["Series"].astype(str)

    subject_list = sorted(pd.unique(df.Participant.ravel()))
    result = {}
    func = getattr(sys.modules[__name__], qap_type)

    # Identify failed subjects
    #failed = ['%s (%s_%s)' % (s['id'], s['session'], s['scan'])
    #          for s in res_dict if 'failed' in s['status']]

    pdf_group = []

    # Generate summary page
    #out_sum = op.join(out_dir, run_name, 'summary_group.pdf')
    #summary_cover(
    #    (qap_type,
    #     datetime.datetime.now().strftime("%Y-%m-%d, %H:%M"),
    #     ", ".join(failed) if len(failed) > 0 else "none"),
    #    is_group=True, out_file=out_sum)
    #pdf_group.append(out_sum)

    # Generate group report
    qc_group = op.join(os.getcwd(), 'qc_measures_group.pdf')
    # Generate violinplots. If successful, add documentation.
    func(df, out_file=qc_group)
    pdf_group.append(qc_group)

    # Generate documentation page
    doc = op.join(os.getcwd(), 'documentation.pdf')

    # Let documentation page fail
    get_documentation(qap_type, doc)
    if doc is not None:
        pdf_group.append(doc)

    if len(pdf_group) > 0:
        out_group_file = op.join(os.getcwd(), '%s_group.pdf' % qap_type)
        # Generate final report with collected pdfs in plots
        concat_pdf(pdf_group, out_group_file)
        result['group'] = {'success': True, 'path': out_group_file}

    if full_reports:
        idx = 1
        # Generate individual reports for subjects
        for subid in [str(sub) for sub in subject_list]:
            print "Generating report for %s.. (%d/%d)" \
                  % (subid, idx, len(subject_list))
            # Get subject-specific info
            subdf = df.loc[df['Participant'] == subid]
            sessions = sorted(pd.unique(subdf.Session.ravel()))
            plots = []
            sess_scans = []
            # Re-build mosaic location
            for sesid in [str(sess) for sess in sessions]:
                sesdf = subdf.loc[subdf['Session'] == sesid]
                scans = sorted(pd.unique(sesdf.Series.ravel()))

                # Each scan has a volume and (optional) fd plot
                for scanid in [str(scan) for scan in scans]:
                    sub_info = [subid, sesid, scanid]
                    sub_path = op.join(out_dir, run_name, '/'.join(sub_info))
                    m = op.join(sub_path, 'qap_mosaic', 'mosaic.pdf')

                    if op.isfile(m):
                        plots.append(m)

                    fd = op.join(sub_path, 'qap_fd', 'fd.pdf')
                    if 'functional_temporal' in qap_type and op.isfile(fd):
                        plots.append(fd)

                sess_scans.append('%s (%s)' % (sesid, ', '.join(scans)))

            #failed = ['%s (%s)' % (s['Session'], s['Series'])
            #          for s in res_dict if 'failed' in s['status'] and
            #          subid in s['id']]

            # Summary cover
            #out_sum = op.join(out_dir, run_name, 'summary_%s.pdf' % subid)
            #summary_cover(
            #    (subid, subid, qap_type,
            #     datetime.datetime.now().strftime("%Y-%m-%d, %H:%M"),
            #     ", ".join(sess_scans),
            #     ",".join(failed) if len(failed) > 0 else "none"),
            #    out_file=out_sum)
            #plots.insert(0, out_sum)

            # Summary (violinplots) of QC measures
            qc_ms = op.join(os.getcwd(), '%s_measures.pdf' % qap_type)
            func(df, subject=subid, out_file=qc_ms)
            plots.append(qc_ms)

            if len(plots) > 0:
                if doc is not None:
                    plots.append(doc)

                # Generate final report with collected pdfs in plots
                sub_path = out_file % subid
                concat_pdf(plots, sub_path)
                result[subid] = {'success': True, 'path': sub_path}

            idx += 1

    return result


def get_documentation(doc_type, out_file):
    import codecs
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


def summary_cover(data, is_group=False, out_file=None):
    import codecs
    from xhtml2pdf import pisa

    # open output file for writing (truncated binary)
    result = open(out_file, "w+b")

    html_file = 'cover_group.html' if is_group else 'cover_subj.html'

    html_dir = op.abspath(
        op.join(op.dirname(__file__), 'html', html_file))

    with codecs.open(html_dir, mode='r', encoding='utf-8') as f:
        html = f.read()

    # convert HTML to PDF
    status = pisa.pisaDocument(html % data, result, encoding='UTF-8')
    result.close()

    # return True on success and False on errors
    return status.err


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
    sessions = sorted(pd.unique(df.Session.ravel()))
    for ss in sessions:
        sesdf = df.copy().loc[df['Session'] == ss]
        scans = pd.unique(sesdf.Series.ravel())
        if sc_split:
            for sc in scans:
                subset = sesdf.loc[sesdf['Series'] == sc]
                if len(subset.index) > 1:
                    if sub_id is None:
                        subtitle = '(%s_%s)' % (ss, sc)
                    else:
                        subtitle = '(Participant %s_%s_%s)' % (sub_id, ss, sc)
                    if condensed:
                        fig = plot_all(sesdf, groups, subject=sub_id,
                                       title='QC measures ' + subtitle)
                    else:
                        fig = plot_measures(
                            sesdf, headers, subject=sub_id,
                            title='QC measures ' + subtitle)
                    if not fig:
                        # this happens if there is a sub_id, but the sub_id
                        # does not have one of the sessions in "sessions"
                        continue
                    report.savefig(fig, dpi=300)
                    fig.clf()
        else:
            if len(sesdf.index) > 1:
                if sub_id is None:
                    subtitle = '(%s)' % (ss)
                else:
                    subtitle = '(Participant %s_%s)' % (sub_id, ss)
                if condensed:
                    fig = plot_all(sesdf, groups, subject=sub_id,
                                   title='QC measures ' + subtitle)
                else:
                    fig = plot_measures(
                        sesdf, headers, subject=sub_id,
                        title='QC measures ' + subtitle)
                if not fig:
                    # this happens if there is a sub_id, but the sub_id does
                    # not have one of the sessions in "sessions"
                    continue
                report.savefig(fig, dpi=300)
                fig.clf()

    report.close()
    plt.close()

    return out_file


def _write_all_reports(df, groups, sc_split=False, condensed=True,
                       out_file='report.pdf'):

    outlist = []
    _write_report(
        df, groups, sc_split=sc_split, condensed=condensed, out_file=out_file)

    subject_list = sorted(pd.unique(df.Participant.ravel()))
    for sub_id in subject_list:
        tpl, _ = op.splitext(op.basename(out_file))
        tpl = op.join(op.dirname(out_file), tpl) + '_%s.pdf'
        outlist.append(_write_report(
            df, groups, sub_id=sub_id, sc_split=sc_split, condensed=condensed,
            out_file=tpl % sub_id))
    return out_file, outlist


def all_anatomical(df, sc_split=False, condensed=True,
                   out_file='anatomical.pdf'):
    groups = [['CNR'],
              ['Cortical Contrast'],
              ['EFC'],
              ['FBER'],
              ['FWHM', 'FWHM_x', 'FWHM_y', 'FWHM_z'],
              ['Qi1'],
              ['SNR']]
    return _write_all_reports(
        df, groups, sc_split=sc_split,
        condensed=condensed, out_file=out_file)


def all_func_temporal(df, sc_split=False, condensed=True,
                      out_file='func_temporal.pdf'):
    groups = [['Fraction of Outliers (Mean)', 'Fraction of Outliers (Median)',
               'Fraction of Outliers (Std Dev)', 'Fraction of Outliers IQR'],
              ['GCOR'],
              ['Quality (Mean)', 'Quality (Median)', 'Quality (Std Dev)',
               'Quality IQR', 'Quality percent outliers'],
              ['RMSD (Mean)', 'RMSD (Median)', 'RMSD (Std Dev)', 'RMSD IQR'],
              ['Std. DVARS (Mean)', 'Std. DVARS (Median)',
               'Std. DVARS percent outliers', 'Std. DVARs IQR']]
    return _write_all_reports(
        df, groups, sc_split=sc_split,
        condensed=condensed, out_file=out_file)


def all_func_spatial(df, sc_split=False, condensed=False,
                     out_file='func_spatial.pdf'):
    groups = [['EFC'],
              ['FBER'],
              ['FWHM', 'FWHM_x', 'FWHM_y', 'FWHM_z'],
              ['Ghost_%s' % a for a in ['x', 'y', 'z']],
              ['SNR']]
    return _write_all_reports(
        df, groups, sc_split=sc_split,
        condensed=condensed, out_file=out_file)


def qap_anatomical_spatial(
        df, subject=None, sc_split=False, condensed=True,
        out_file='anatomical.pdf'):
    groups = [['CNR'],
              ['Cortical Contrast'],
              ['EFC'],
              ['FBER'],
              ['FWHM', 'FWHM_x', 'FWHM_y', 'FWHM_z'],
              ['Qi1'],
              ['SNR']]
    return _write_report(
        df, groups, sub_id=subject, sc_split=sc_split, condensed=condensed,
        out_file=out_file)


def qap_functional_temporal(
        df, subject=None, sc_split=False, condensed=True,
        out_file='func_temporal.pdf'):
    groups = [['Fraction of Outliers (Mean)', 'Fraction of Outliers (Median)',
               'Fraction of Outliers (Std Dev)', 'Fraction of Outliers IQR'],
              ['GCOR'],
              ['Quality (Mean)', 'Quality (Median)', 'Quality (Std Dev)',
               'Quality IQR', 'Quality percent outliers'],
              ['RMSD (Mean)', 'RMSD (Median)', 'RMSD (Std Dev)', 'RMSD IQR'],
              ['Std. DVARS (Mean)', 'Std. DVARS (Median)',
               'Std. DVARS percent outliers', 'Std. DVARs IQR']]
    return _write_report(
        df, groups, sub_id=subject, sc_split=sc_split, condensed=condensed,
        out_file=out_file)


def qap_functional_spatial(
        df, subject=None, sc_split=False, condensed=True,
        out_file='func_spatial.pdf'):
    groups = [['EFC'],
              ['FBER'],
              ['FWHM', 'FWHM_x', 'FWHM_y', 'FWHM_z'],
              ['Ghost_%s' % a for a in ['x', 'y', 'z']],
              ['SNR']]
    return _write_report(
        df, groups, sub_id=subject, sc_split=sc_split, condensed=condensed,
        out_file=out_file)
