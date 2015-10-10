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

from .plotting import (plot_mosaic, plot_dist)

from nipype import logging
iflogger = logging.getLogger('interface')


class PlotMosaicInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True,
                   desc='File to be plotted')
    in_mask = File(exists=True, desc='Overlay mask')
    title = traits.Str('Volume', desc='modality name to be prepended')
    subject = traits.Str(desc='Subject id')
    metadata = traits.List(traits.Str, desc='additional metadata')
    figsize = traits.Tuple(traits.Float, traits.Float, default=(11.69, 8.27),
                           desc='Figure size')
    dpi = traits.Int(300, desc='Desired DPI of figure')
    out_file = File('mosaic.pdf', desc='output file name')


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

        fig = plot_mosaic(
            self.inputs.in_file,
            title=title,
            overlay_mask=mask,
            figsize=self.inputs.figsize)

        fig.savefig(self.inputs.out_file, dpi=self.inputs.dpi)

        return runtime

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = op.abspath(self.inputs.out_file)
        return outputs
