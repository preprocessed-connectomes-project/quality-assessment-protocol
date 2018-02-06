#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import os
import os.path as op

import nibabel as nb
import numpy as np

from nipype.interfaces.base import (BaseInterface, traits, TraitedSpec, File,
                                    InputMultiPath, OutputMultiPath,
                                    BaseInterfaceInputSpec, isdefined,
                                    DynamicTraitedSpec, Undefined)

from .plotting import (plot_mosaic, grayplot)

from nipype import logging
iflogger = logging.getLogger('interface')


class PlotMosaicInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True,
                   desc='File to be plotted')
    in_mask = File(exists=True, desc='Overlay mask')
    title = traits.Str('Volume', usedefault=True,
                       desc='modality name to be prepended')
    subject = traits.Str(desc='Subject id')
    metadata = traits.List(traits.Str, desc='additional metadata')
    figsize = traits.Tuple(
        (11.69, 8.27), traits.Float, traits.Float, usedefault=True,
        desc='Figure size')
    dpi = traits.Int(300, usedefault=True, desc='Desired DPI of figure')
    out_file = File('mosaic.png', usedefault=True, desc='output file name')


class PlotMosaicOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='output pdf file')


class PlotMosaic(BaseInterface):

    """
    Plots slices of a 3D volume into a pdf file
    """
    input_spec = PlotMosaicInputSpec
    output_spec = PlotMosaicOutputSpec

    def _run_interface(self, runtime):
        mask = None
        if isdefined(self.inputs.in_mask):
            mask = self.inputs.in_mask

        title = self.inputs.title
        if isdefined(self.inputs.subject):
            title += ', subject %s' % self.inputs.subject

        if isdefined(self.inputs.metadata):
            title += ' (' + '_'.join(self.inputs.metadata) + ')'

        if isdefined(self.inputs.figsize):
            fig = plot_mosaic(
                self.inputs.in_file,
                title=title,
                overlay_mask=mask,
                figsize=self.inputs.figsize)
        else:
            fig = plot_mosaic(
                self.inputs.in_file,
                title=title,
                overlay_mask=mask)

        fig.savefig(self.inputs.out_file, dpi=self.inputs.dpi)

        return runtime

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = op.abspath(self.inputs.out_file)
        return outputs


class PlotFDInputSpec(BaseInterfaceInputSpec):
    meanfd_file = File(exists=True, mandatory=True, desc='mean fd file to be plotted')
    dvars = traits.Dict(mandatory=True, desc='dvars float array be plotted')
    global_signal = traits.List(traits.Float, mandatory=True, desc='global signal to be plotted')
    metadata = traits.List(traits.Str, mandatory=True, desc='additional metadata')
    title = traits.Str('Mean FD, DVARS ad global Signal', usedefault=True,
                       desc='modality name to be prepended')
    subject = traits.Str(desc='Subject id')
    dpi = traits.Int(300, usedefault=True, desc='Desired DPI of figure')
    out_file = File('fd.png', usedefault=True, desc='output file name')


class PlotFDOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='output pdf file')


class PlotFD(BaseInterface):

    """
    Plots the frame displacement of a dataset
    """
    input_spec = PlotFDInputSpec
    output_spec = PlotFDOutputSpec

    def _run_interface(self, runtime):
        title = self.inputs.title
        if isdefined(self.inputs.subject):
            title += ', subject %s' % self.inputs.subject

        fig = plot_fd(self.inputs.meanfd_file, self.inputs.dvars,
                      self.inputs.global_signal, self.inputs.metadata)

        fig.savefig(self.inputs.out_file, dpi=float(self.inputs.dpi))

        return runtime

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = op.abspath(self.inputs.out_file)
        return outputs

class GrayPlotInputSpec(BaseInterfaceInputSpec):
    func_file = File(exists=True, mandatory=True, desc='functional file to be plotted in gray plot')
    mask_file = File(exists=True, mandatory=True, desc='file for mask functional data')
    meanfd_file = File(exists=True, mandatory=True, desc='mean fd file to be plotted')
    dvars = traits.List(traits.Float, mandatory=True, desc='dvars float array be plotted')
    global_signal = traits.List(traits.Float, mandatory=True, desc='global signal to be plotted')
    metadata = traits.List(traits.Str, mandatory=True, desc='additional metadata')
    title = traits.Str('Timeseries  Plot', usedefault=True,
                       desc='modality name to be prepended')
    subject = traits.Str(desc='Subject id')
    dpi = traits.Int(300, usedefault=True, desc='Desired DPI of figure')
    out_file = File('grayplot.png', usedefault=True, desc='output file name')
    out_cluster = File('grayplot_cluster.nii.gz', usedefault=True, 
        desc='output file name of the cluster generated using grayplot')


class GrayPlotOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='output png file')
    out_cluster = File(exists=True, desc='cluster generated using grayplot')


class GrayPlot(BaseInterface):
    """
    Plots the frame displacement of a dataset
    """
    input_spec = GrayPlotInputSpec
    output_spec = GrayPlotOutputSpec

    def _run_interface(self, runtime):
        title = self.inputs.title
        if isdefined(self.inputs.subject):
            title += ', %s' % self.inputs.subject
        fig, cluser = grayplot(self.inputs.func_file, self.inputs.mask_file, self.inputs.meanfd_file, self.inputs.dvars, self.inputs.global_signal, self.inputs.metadata, title=title)

        fig.savefig(self.inputs.out_file, dpi=float(self.inputs.dpi))
        nb.save(cluser, self.inputs.out_cluster)

        return runtime

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = op.abspath(self.inputs.out_file)
        return outputs
