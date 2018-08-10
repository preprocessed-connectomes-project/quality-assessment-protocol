#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import os
import copy

import nipype.pipeline.engine as nipype_pipe_engine
import nipype.interfaces.utility as nipype_utility
from nipype.interfaces.afni import utils as nipype_afni_utils
from nipype.interfaces.afni import preprocess as nipype_afni_preprocess

import anatomical_preproc
import functional_preproc
import qap_utils
import qap_workflows_utils
import spatial_qc
import temporal_qc
import qap.viz.interfaces as qap_viz


def qap_mask_workflow(workflow, resource_pool, config, name="_"):
    """Build and run a Nipype workflow to create the QAP anatomical head mask.

    - If any resources/outputs required by this workflow are not in the
      resource pool, this workflow will call pre-requisite workflow builder
      functions to further populate the pipeline with workflows which will
      calculate/generate these necessary pre-requisites.

    Expected Resources in Resource Pool:
      - anatomical_reorient: The deobliqued, reoriented anatomical scan.
      - allineate_linear_xfm: The linear registration transform matrix from
                              AFNI's 3dAllineate.

    New Resources Added to Resource Pool
      - qap_head_mask: A binary mask of the head and the region in front of
                       the mouth and below the nose.
      - whole_head_mask: A binary mask of only the head.
      - skull_only_mask: A binary mask of the head minus the slice mask.

    Workflow Steps
      1. AFNI 3dClipLevel on anatomical reorient to get a threshold value.
      2. create_expr_string function node to generate the arg for AFNI 3dcalc.
      3. AFNI 3dcalc to create a binary mask of anatomical reorient using the
         threshold value from 3dClipLevel.
      4. AFNI 3dmask_tool to dilate and erode the binary mask six times each
         way to remove gaps and holes.
      5. slice_head_mask function node to create the binary mask for the nose/
         mouth region.
      6. 3dcalc to add the slice mask to the head mask, and to also subtract
         the slice mask from the head mask, to create the "skull_only_mask".

    :type workflow: Nipype workflow object
    :param workflow: A Nipype workflow object which can already contain other
                     connected nodes; this function will insert the following
                     workflow into this one provided.
    :type resource_pool: dict
    :param resource_pool: A dictionary defining input files and pointers to
                          Nipype node outputs / workflow connections; the keys
                          are the resource names.
    :type config: dict
    :param config: A dictionary defining the configuration settings for the
                   workflow, such as directory paths or toggled options.
    :type name: str
    :param name: (default: "_") A string to append to the end of each node
                 name.
    :rtype: Nipype workflow object
    :return: The Nipype workflow originally provided, but with this function's
              sub-workflow connected into it.
    :rtype: dict
    :return: The resource pool originally provided, but updated (if
             applicable) with the newest outputs and connections.
    """

    if 'anat_linear_xfm' not in resource_pool.keys():
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = anatomical_preproc.afni_anatomical_linear_registration(workflow, resource_pool,
                                                                                         config, name)

        if resource_pool == old_rp:
            return workflow, resource_pool

    if 'anat_reorient' not in resource_pool.keys():
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            anatomical_preproc.anatomical_reorient_workflow(workflow, resource_pool, config, name)

        if resource_pool == old_rp:
            return workflow, resource_pool

    # find the clipping level for thresholding the head mask
    clip_level = nipype_pipe_engine.Node(interface=nipype_afni_preprocess.ClipLevel(),
                                         name='qap_headmask_clip_level{}'.format(name))

    # create the expression string for Calc in mask_skull node
    create_expr_string = nipype_pipe_engine.Node(nipype_utility.Function(
        input_names=['clip_level_value'],
        output_names=['expr_string'],
        function=qap_utils.create_expr_string),
        name='qap_headmask_create_expr_string{0}'.format(name))

    workflow.connect(clip_level, 'clip_val', create_expr_string, 'clip_level_value')

    # let's create a binary mask of the skull image with that threshold
    try:
        mask_skull = nipype_pipe_engine.Node(interface=nipype_afni_utils.Calc(),
                                             name='qap_headmask_mask_skull{0}'.format(name))
    except AttributeError:
        mask_skull = nipype_pipe_engine.Node(interface=nipype_afni_utils.Calc(),
                                             name='qap_headmask_mask_skull{0}'.format(name))

    mask_skull.inputs.outputtype = "NIFTI_GZ"

    workflow.connect(create_expr_string, 'expr_string', mask_skull, 'expr')

    try:
        dilate_erode = nipype_pipe_engine.Node(interface=nipype_afni_utils.MaskTool(),
                                               name='qap_headmask_mask_tool{0}'.format(name))
    except AttributeError:
        dilate_erode = nipype_pipe_engine.Node(interface=nipype_afni_utils.MaskTool(),
                                               name='qap_headmask_mask_tool{0}'.format(name))

    dilate_erode.inputs.dilate_inputs = "6 -6"
    dilate_erode.inputs.outputtype = "NIFTI_GZ"

    workflow.connect(mask_skull, 'out_file', dilate_erode, 'in_file')

    # create Scipy hole-filling function node
    fill = nipype_pipe_engine.Node(nipype_utility.Function(input_names=['in_file'],
                                                           output_names=['out_file'],
                                                           function=qap_workflows_utils.scipy_fill),
                                   name='qap_headmask_scipy_fill{0}'.format(name))

    workflow.connect(dilate_erode, 'out_file', fill, 'in_file')

    slice_head_mask = nipype_pipe_engine.Node(nipype_utility.Function(
        input_names=['infile', 'transform', 'standard'],
        output_names=['outfile_path'], function=qap_workflows_utils.slice_head_mask),
        name='qap_headmask_slice_head_mask{0}'.format(name))

    try:
        combine_masks = nipype_pipe_engine.Node(interface=nipype_afni_utils.Calc(),
                                                name='qap_headmask_combine_masks{0}'.format(name))
    except AttributeError:
        combine_masks = nipype_pipe_engine.Node(interface=nipype_afni_utils.Calc(),
                                                name='qap_headmask_combine_masks{0}'.format(name))

    combine_masks.inputs.expr = "(a+b)-(a*b)"
    combine_masks.inputs.outputtype = "NIFTI_GZ"

    # subtract the slice mask from the original head mask to create a
    # skull-only mask for FG calculations
    try:
        subtract_mask = nipype_pipe_engine.Node(interface=nipype_afni_utils.Calc(),
                                                name='qap_headmask_subtract_masks{0}'.format(name))
    except AttributeError:
        subtract_mask = nipype_pipe_engine.Node(interface=nipype_afni_utils.Calc(),
                                                name='qap_headmask_subtract_masks{0}'.format(name))

    subtract_mask.inputs.expr = "a-b"
    subtract_mask.inputs.outputtype = "NIFTI_GZ"

    if isinstance(resource_pool['anat_reorient'], tuple):
        node, out_file = resource_pool['anat_reorient']
        workflow.connect([
            (node, clip_level,      [(out_file, 'in_file')]),
            (node, mask_skull,      [(out_file, 'in_file_a')]),
            (node, slice_head_mask, [(out_file, 'infile')])
        ])
    else:
        clip_level.inputs.in_file = resource_pool['anat_reorient']
        mask_skull.inputs.in_file_a = resource_pool['anat_reorient']
        slice_head_mask.inputs.infile = resource_pool['anat_reorient']

    if isinstance(resource_pool['anat_linear_xfm'], tuple):
        node, out_file = resource_pool['anat_linear_xfm']
        workflow.connect(node, out_file, slice_head_mask, 'transform')
    else:
        slice_head_mask.inputs.transform = \
            resource_pool['anat_linear_xfm']

    workflow.connect([
        (fill, combine_masks, [('out_file', 'in_file_a')]),
        (slice_head_mask, combine_masks, [('outfile_path', 'in_file_b')])
    ])

    workflow.connect(combine_masks, 'out_file', subtract_mask, 'in_file_a')
    workflow.connect(slice_head_mask, 'outfile_path',
                     subtract_mask, 'in_file_b')

    resource_pool['anat_qap_head_mask'] = (combine_masks, 'out_file')
    resource_pool['anat_whole_head_mask'] = (fill, 'out_file')
    resource_pool['anat_skull_only_mask'] = (subtract_mask, 'out_file')

    return workflow, resource_pool


def qap_gather_header_info(workflow, resource_pool, config, name="_", data_type="anat"):
    """Build and run a Nipype workflow to extract the NIFTI header information
    from an input file and insert it into a dictionary.

    - If any resources/outputs required by this workflow are not in the
      resource pool, this workflow will call pre-requisite workflow builder
      functions to further populate the pipeline with workflows which will
      calculate/generate these necessary pre-requisites.

    Expected Resources in Resource Pool
      - anatomical_scan: (optional) The raw anatomical scan.
      - functional_scan: (optional) The raw functional scan.

    New Resources Added to Resource Pool
      - anatomical_header_info: (if anatomical scan provided) Dictionary of
                                NIFTI file's header information.
      - functional_header_info: (if functional scan provided) Dictionary of
                                NIFTI file's header information.

    Workflow Steps
      1. create_header_dict_entry function node to generate the header info
         dictionary.
      2. write_json function node to write the information to the output JSON
         file (or update an already existing output JSON file).

    :type workflow: Nipype workflow object
    :param workflow: A Nipype workflow object which can already contain other
                     connected nodes; this function will insert the following
                     workflow into this one provided.
    :type resource_pool: dict
    :param resource_pool: A dictionary defining input files and pointers to
                          Nipype node outputs / workflow connections; the keys
                          are the resource names.
    :type config: dict
    :param config: A dictionary defining the configuration settings for the
                   workflow, such as directory paths or toggled options.
    :type name: str
    :param name: (default: "_") A string to append to the end of each node
                 name.

    :type data_type: str
    :param data_type: The BIDS type for the data being extracted, either 'anat' or 'func'

    :rtype: Nipype workflow object
    :return: The Nipype workflow originally provided, but with this function's
              sub-workflow connected into it.
    :rtype: dict
    :return: The resource pool originally provided, but updated (if
             applicable) with the newest outputs and connections.
    """

    gather_imports = ['import os',
                      'import nibabel as nb']

    gather_header = nipype_pipe_engine.Node(
        nipype_utility.Function(
            input_names=['in_file', 'subject', 'session', 'scan', 'type'],
            output_names=['qap_dict'],
            function=qap_workflows_utils.create_header_dict_entry,
            imports=gather_imports
        ),
        name="gather_header_info_{0}{1}".format(data_type, name)
    )
    
    gather_header.inputs.subject = config["subject_id"]
    gather_header.inputs.session = config["session_id"]
    gather_header.inputs.scan = config["scan_id"]

    if data_type == "anat":
        if "anatomical_scan" in resource_pool.keys():
            gather_header.inputs.in_file = resource_pool["anatomical_scan"]
            gather_header.inputs.type = data_type
        else:
            raise ValueError("Could not find anatomical_scan in resource pool for data type {0}.".format(data_type))

    elif "func" in data_type:
        if "functional_scan" in resource_pool.keys():
            gather_header.inputs.in_file = resource_pool["functional_scan"]
            gather_header.inputs.type = data_type
        else:
            raise ValueError("Could not find anatomical_scan in resource pool for data type {0}.".format(data_type))

    else:
        raise ValueError("Do not know how to handle data type {0}".format(data_type))

    out_dir = os.path.join(config['output_directory'], config["run_name"],
                           config["site_name"], config["subject_id"],
                           config["session_id"])

    out_json = os.path.join(out_dir, "{0}_{1}_{2}_qap-{3}.json".format(config["subject_id"], config["session_id"],
                                                                       config["scan_id"], data_type))
    header_to_json = nipype_pipe_engine.Node(nipype_utility.Function(
                                 input_names=["output_dict",
                                              "json_file"],
                                 output_names=["json_file"],
                                 function=qap_utils.write_json),
                             name="qap_header_to_json_{0}{1}".format(data_type, name))
    header_to_json.inputs.json_file = out_json

    workflow.connect(gather_header, 'qap_dict', header_to_json, 'output_dict')
    resource_pool['{0}_header_info'.format(data_type)] = out_json

    return workflow, resource_pool


def calculate_artifacts_background(workflow, resource_pool, config, name="_"):
    """

    Creates workflow for calculating the fraction of artifact voxels and inserts it into the pipeline.

    :param workflow: the workflow to which the newly created workflow will be appended
    :param resource_pool: the resource pool corresponding to the workflow
    :param config: pipeline configuration dictionary
    :param name: name assigned to the created workflow
    :return: tuple (workflow, resource_pool) after the new workflow has been incorporated
    """

    if 'anat_qap_head_mask' not in resource_pool.keys():
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            qap_mask_workflow(workflow, resource_pool, config, name)
        if resource_pool == old_rp:
            return workflow, resource_pool

    if 'anat_reorient' not in resource_pool.keys():
        old_rp = copy.copy(resource_pool)
        workflow, new_resource_pool = \
            anatomical_preproc.anatomical_reorient_workflow(workflow, resource_pool, config,
                                                            name)
        if resource_pool == old_rp:
            return workflow, resource_pool

    calculate_artifacts = \
        nipype_pipe_engine.Node(nipype_utility.Function(input_names=['anatomical_reorient',
                                                                     'qap_head_mask_path',
                                                                     'exclude_zeroes'],
                                                        output_names=['fav_bg_file', 'bg_mask_file'],
                                                        function=spatial_qc.artifacts),
                                name="calculate_artifacts{0}".format(name))

    calculate_artifacts.inputs.exclude_zeroes = config["exclude_zeros"]

    if isinstance(resource_pool["anat_reorient"], tuple):
        node, out_file = resource_pool["anat_reorient"]
        workflow.connect(node, out_file,
                         calculate_artifacts, 'anatomical_reorient')
    else:
        calculate_artifacts.inputs.anatomical_reorient = \
            resource_pool["anat_reorient"]

    if isinstance(resource_pool["anat_qap_head_mask"], tuple):
        node, out_file = resource_pool["anat_qap_head_mask"]
        workflow.connect(node, out_file,
                         calculate_artifacts, 'qap_head_mask_path')
    else:
        calculate_artifacts.inputs.qap_head_mask_path = \
            resource_pool["anat_qap_head_mask"]

    resource_pool["anat_fav_artifacts_background"] = \
        (calculate_artifacts, 'fav_bg_file')
    resource_pool["anat_qap_bg_head_mask"] = \
        (calculate_artifacts, 'bg_mask_file')

    return workflow, resource_pool


def calculate_temporal_std(workflow, resource_pool, config=None, name="_"):
    """
    Creates a workflow the calculate temporal standard deviation image and inserts into the passed in workflow and
    resource pool

    :param workflow: workflow to which the newly created workflow will be appended
    :param resource_pool: resource pool to which the workflow outputs will be added
    :param config: pipeline configuration dictionary, not used for anything but included to make consistent with other
        workflow creation functions
    :param name: name given to the created workflow
    :return: tuple (workflow, resource_pool) with the new workflow added
    """

    calculate_tstd_map = \
        nipype_pipe_engine.Node(nipype_utility.Function(input_names=['func_reorient',
                                                                     'func_mask'],
                                                        output_names=['temporal_std_map_file'],
                                                        function=qap_workflows_utils.get_temporal_std_map),
                                name="calculate_temporal_std{0}".format(name))

    if isinstance(resource_pool["func_reorient"], tuple):
        node, out_file = resource_pool["func_reorient"]
        workflow.connect(node, out_file, calculate_tstd_map, 'func_reorient')
    else:
        calculate_tstd_map.inputs.func_reorient = \
            resource_pool["func_reorient"]

    if isinstance(resource_pool["func_brain_mask"], tuple):
        node, out_file = resource_pool["func_brain_mask"]
        workflow.connect(node, out_file, calculate_tstd_map, 'func_mask')
    else:
        calculate_tstd_map.inputs.func_mask = \
            resource_pool["func_brain_mask"]

    resource_pool['func_temporal_std_map'] = (calculate_tstd_map,
                                              'temporal_std_map_file')

    return workflow, resource_pool


def calculate_sfs_workflow(workflow, resource_pool, config, name="_"):
    """ create workflow for calculating signal fluctuation sensitivity and insert into passed in workflow and resource
        pool

    :param workflow: workflow to which the newly created workflow will be appended
    :param resource_pool: resources pool to which the outputs of the new workflow will be added
    :param config: pipeline configuration dictionary
    :param name: name given to the created workflow
    :return: tuple (workflow, resource_pool) with the new workflow included
    """

    if 'func_mean' not in resource_pool.keys():
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            functional_preproc.mean_functional_workflow(workflow, resource_pool, config, name)
        if resource_pool == old_rp:
            return workflow, resource_pool
            
    if 'func_temporal_std_map' not in resource_pool.keys():
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            calculate_temporal_std(workflow, resource_pool, config, name)
        if resource_pool == old_rp:
            return workflow, resource_pool

    calculate_sfs = nipype_pipe_engine.Node(nipype_utility.Function(input_names=['func',
                                                                                 'func_mean',
                                                                                 'func_mask',
                                                                                 'temporal_std_file'],
                                                                    output_names=['sfs_file',
                                                                                  'est_nuisance_file'],
                                                                    function=qap_workflows_utils.sfs_timeseries),
                                            name="calculate_sfs{0}".format(name))

    if isinstance(resource_pool["func_reorient"], tuple):
        node, out_file = resource_pool["func_reorient"]
        workflow.connect(node, out_file, calculate_sfs, 'func')
    else:
        calculate_sfs.inputs.func = \
            resource_pool["func_reorient"]

    if isinstance(resource_pool["func_mean"], tuple):
        node, out_file = resource_pool["func_mean"]
        workflow.connect(node, out_file, calculate_sfs, 'func_mean')
    else:
        calculate_sfs.inputs.func_mean = resource_pool["func_mean"]

    if isinstance(resource_pool["func_brain_mask"], tuple):
        node, out_file = resource_pool["func_brain_mask"]
        workflow.connect(node, out_file, calculate_sfs, 'func_mask')
    else:
        calculate_sfs.inputs.func_mask = \
            resource_pool["func_brain_mask"]

    if isinstance(resource_pool["func_temporal_std_map"], tuple):
        node, out_file = resource_pool["func_temporal_std_map"]
        workflow.connect(node, out_file, calculate_sfs, 'temporal_std_file')
    else:
        calculate_sfs.inputs.temporal_std_file = \
            resource_pool["func_temporal_std_map"]

    resource_pool['func_SFS'] = (calculate_sfs, 'sfs_file')
    resource_pool['func_estimated_nuisance'] = \
        (calculate_sfs, 'est_nuisance_file')

    return workflow, resource_pool


def qap_anatomical_spatial_workflow(workflow, resource_pool, config, name="_"):
    """Build and run a Nipype workflow to calculate the QAP anatomical spatial
    quality measures.

    - If any resources/outputs required by this workflow are not in the
      resource pool, this workflow will call pre-requisite workflow builder
      functions to further populate the pipeline with workflows which will
      calculate/generate these necessary pre-requisites.

    Expected Resources in Resource Pool
      - anatomical_reorient: The deobliqued, reoriented anatomical scan.
      - qap_head_mask: A binary mask of the head and the region in front of
                       the mouth and below the nose.
      - anatomical_gm_mask: The binary tissue segmentation map for gray
                            matter.
      - anatomical_wm_mask: The binary tissue segmentation map for white
                            matter.
      - anatomical_csf_mask: The binary tissue segmentation map for CSF.

    New Resources Added to Resource Pool
      - qap_anatomical_spatial: The path to the output JSON file containing
                                the participant's QAP measure values.
      - anat_spat_csv: The path to the CSV file containing the QAP measure
                       values.
      - qap_mosaic: (if enabled) The path to the mosaic QC report file.

    Workflow Steps:
      1. qap_anatomical_spatial function node to calculate the QAP measures.
      2. PlotMosaic() node (if enabled) to generate QC mosaic.
      3. qap_anatomical_spatial_to_json function node to write/update numbers
         to the output JSON file.

    :type workflow: Nipype workflow object
    :param workflow: A Nipype workflow object which can already contain other
                     connected nodes; this function will insert the following
                     workflow into this one provided.
    :type resource_pool: dict
    :param resource_pool: A dictionary defining input files and pointers to
                          Nipype node outputs / workflow connections; the keys
                          are the resource names.
    :type config: dict
    :param config: A dictionary defining the configuration settings for the
                   workflow, such as directory paths or toggled options.
    :type name: str
    :param name: (default: "_") A string to append to the end of each node
                 name.

    :rtype: Nipype workflow object
    :return: The Nipype workflow originally provided, but with this function's
              sub-workflow connected into it.
    :rtype: dict
    :return: The resource pool originally provided, but updated (if
             applicable) with the newest outputs and connections.
    """

    qap_utils.check_config_settings(config, "anatomical_template")

    if "exclude_zeros" not in config.keys():
        config["exclude_zeros"] = False

    if 'anat_fav_artifacts_background' not in resource_pool.keys():
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            calculate_artifacts_background(workflow, resource_pool, config,
                                           name)
        if resource_pool == old_rp:
            return workflow, resource_pool

    if ('anat_gm_mask' not in resource_pool.keys()) or \
            ('anat_wm_mask' not in resource_pool.keys()) or \
            ('anat_csf_mask' not in resource_pool.keys()):

        old_rp = copy.copy(resource_pool)
        workflow, new_resource_pool = \
            anatomical_preproc.afni_segmentation_workflow(workflow, resource_pool, config, name)

        if resource_pool == old_rp:
            return workflow, resource_pool

    spatial = nipype_pipe_engine.Node(nipype_utility.Function(
        input_names=['anatomical_reorient', 'qap_head_mask_path',
                     'qap_bg_head_mask_path', 'whole_head_mask_path',
                     'skull_mask_path', 'anatomical_gm_mask',
                     'anatomical_wm_mask', 'anatomical_csf_mask',
                     'fav_artifacts', 'subject_id', 'session_id', 'scan_id',
                     'run_name', 'site_name', 'exclude_zeroes', 'out_vox',
                     'session_output_dir', 'starter'],
        output_names=['qap'], function=qap_workflows_utils.qap_anatomical_spatial),
        name='qap_anatomical_spatial{0}'.format(name))

    # Subject infos
    spatial.inputs.subject_id = config['subject_id']
    spatial.inputs.session_id = config['session_id']
    spatial.inputs.scan_id = config['scan_id']
    spatial.inputs.run_name = config['run_name']
    spatial.inputs.exclude_zeroes = config['exclude_zeros']
    spatial.inputs.out_vox = True
    spatial.inputs.session_output_dir = config['output_directory']

    node, out_file = resource_pool['starter']
    workflow.connect(node, out_file, spatial, 'starter')

    if 'site_name' in resource_pool.keys():
        spatial.inputs.site_name = resource_pool['site_name']
    elif 'site_name' in config.keys():
        spatial.inputs.site_name = config['site_name']

    if isinstance(resource_pool['anat_reorient'], tuple):
        node, out_file = resource_pool['anat_reorient']
        workflow.connect(node, out_file, spatial, 'anatomical_reorient')
    else:
        spatial.inputs.anatomical_reorient = \
            resource_pool['anat_reorient']

    if isinstance(resource_pool['anat_qap_head_mask'], tuple):
        node, out_file = resource_pool['anat_qap_head_mask']
        workflow.connect(node, out_file, spatial, 'qap_head_mask_path')
    else:
        spatial.inputs.qap_head_mask_path = \
            resource_pool['anat_qap_head_mask']

    if isinstance(resource_pool['anat_qap_bg_head_mask'], tuple):
        node, out_file = resource_pool['anat_qap_bg_head_mask']
        workflow.connect(node, out_file, spatial, 'qap_bg_head_mask_path')
    else:
        spatial.inputs.qap_bg_head_mask_path = \
            resource_pool['anat_qap_bg_head_mask']

    if isinstance(resource_pool['anat_whole_head_mask'], tuple):
        node, out_file = resource_pool['anat_whole_head_mask']
        workflow.connect(node, out_file, spatial, 'whole_head_mask_path')
    else:
        spatial.inputs.whole_head_mask_path = \
            resource_pool['anat_whole_head_mask']

    if isinstance(resource_pool['anat_skull_only_mask'], tuple):
        node, out_file = resource_pool['anat_skull_only_mask']
        workflow.connect(node, out_file, spatial, 'skull_mask_path')
    else:
        spatial.inputs.skull_mask_path = resource_pool['anat_skull_only_mask']

    if isinstance(resource_pool['anat_gm_mask'], tuple):
        node, out_file = resource_pool['anat_gm_mask']
        workflow.connect(node, out_file, spatial, 'anatomical_gm_mask')
    else:
        spatial.inputs.anatomical_gm_mask = \
            resource_pool['anat_gm_mask']

    if len(resource_pool['anat_wm_mask']) == 2:
        node, out_file = resource_pool['anat_wm_mask']
        workflow.connect(node, out_file, spatial, 'anatomical_wm_mask')
    else:
        spatial.inputs.anatomical_wm_mask = \
            resource_pool['anat_wm_mask']

    if len(resource_pool['anat_csf_mask']) == 2:
        node, out_file = resource_pool['anat_csf_mask']
        workflow.connect(node, out_file, spatial, 'anatomical_csf_mask')
    else:
        spatial.inputs.anatomical_csf_mask = \
            resource_pool['anat_csf_mask']

    if isinstance(resource_pool['anat_fav_artifacts_background'], tuple):
        node, out_file = resource_pool['anat_fav_artifacts_background']
        workflow.connect(node, out_file, spatial, 'fav_artifacts')
    else:
        spatial.inputs.fav_artifacts = \
            resource_pool['anat_fav_artifacts_background']

    if config.get('write_report', False):
        plot = nipype_pipe_engine.Node(qap_viz.PlotMosaic(), name='anat_plot_mosaic{0}'.format(name))
        plot.inputs.subject = config['subject_id']

        metadata = [config['session_id'], config['scan_id']]
        if 'site_name' in config.keys():
            metadata.append(config['site_name'])

        plot.inputs.metadata = metadata
        plot.inputs.title = 'Anatomical reoriented'

        if len(resource_pool['anat_reorient']) == 2:
            node, out_file = resource_pool['anat_reorient']
            workflow.connect(node, out_file, plot, 'in_file')
        else:
            plot.inputs.in_file = resource_pool['anat_reorient']

        resource_pool['mean_epi_mosaic'] = (plot, 'out_file')
        resource_pool['qap_mosaic'] = (plot, 'out_file')

    out_dir = os.path.join(config['output_directory'], config["run_name"],
                           config["site_name"], config["subject_id"],
                           config["session_id"])
    out_json = os.path.join(out_dir, "{0}_{1}_{2}_qap-anat.json".format(config["subject_id"], config["session_id"],
                                                                        config["scan_id"]))

    spatial_to_json = nipype_pipe_engine.Node(nipype_utility.Function(
                                  input_names=["output_dict",
                                               "json_file"],
                                  output_names=["json_file"],
                                  function=qap_utils.write_json),
                              name="qap_anatomical_spatial_to_json{0}".format(name))
    spatial_to_json.inputs.json_file = out_json

    workflow.connect(spatial, 'qap', spatial_to_json, 'output_dict')
    resource_pool['qap_anatomical_spatial'] = out_json

    return workflow, resource_pool


def qap_functional_workflow(workflow, resource_pool, config, name="_"):
    """Build and run a Nipype workflow to calculate the QAP functional 
    temporal quality measures.

    - If any resources/outputs required by this workflow are not in the
      resource pool, this workflow will call pre-requisite workflow builder
      functions to further populate the pipeline with workflows which will
      calculate/generate these necessary pre-requisites.

    Expected Resources in Resource Pool
      - func_reorient: The deobliqued, reoriented 4D functional timeseries.
      - functional_brain_mask: A binary mask of the brain in the functional
                               image.
      - inverted_functional_brain_mask: A binary mask of the inversion of the
                                        functional brain mask.
      - coordinate_transformation: The matrix transformation from AFNI's
                                   3dvolreg (--1Dmatrix_save option).
      - mcflirt_rel_rms: (if no coordinate_transformation) The matrix
                         transformation from FSL's Mcflirt.

    New Resources Added to Resource Pool
      - qap_functional_temporal: The path to the output JSON file containing
                                 the participant's QAP measure values.
      - func_temp_csv: The path to the CSV file containing the QAP measure
                       values.
      - qap_mosaic: (if enabled) The path to the mosaic QC report file.

    Workflow Steps
      1. fd_jenkinson function node to calculate RMSD.
      2. qap_functional_temporal function node to calculate the QAP measures.
      3. PlotMosaic() node (if enabled) to generate QC mosaic.
      4. qap_functional_temporal_to_json function node to write/update numbers
         to the output JSON file.

    :type workflow: Nipype workflow object
    :param workflow: A Nipype workflow object which can already contain other
                     connected nodes; this function will insert the following
                     workflow into this one provided.
    :type resource_pool: dict
    :param resource_pool: A dictionary defining input files and pointers to
                          Nipype node outputs / workflow connections; the keys
                          are the resource names.
    :type config: dict
    :param config: A dictionary defining the configuration settings for the
                   workflow, such as directory paths or toggled options.
    :type name: str
    :param name: (default: "_") A string to append to the end of each node
                 name.
    :rtype: Nipype workflow object
    :return: The Nipype workflow originally provided, but with this function's
              sub-workflow connected into it.
    :rtype: dict
    :return: The resource pool originally provided, but updated (if
             applicable) with the newest outputs and connections.
    """

    # ensures functional_brain_mask is created as well
    if 'func_inverted_brain_mask' not in resource_pool.keys():
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            functional_preproc.invert_functional_brain_mask_workflow(workflow, resource_pool,
                                                                     config, name)
        if resource_pool == old_rp:
            return workflow, resource_pool

    if ('func_motion_correct' not in resource_pool.keys()) or \
        ('func_coordinate_transformation' not in resource_pool.keys() and
            'mcflirt_rel_rms' not in resource_pool.keys()):
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            functional_preproc.func_motion_correct_workflow(workflow, resource_pool, config,
                                                            name)
        if resource_pool == old_rp:
            return workflow, resource_pool

    if 'func_mean' not in resource_pool.keys():
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            functional_preproc.mean_functional_workflow(workflow, resource_pool, config, name)
        if resource_pool == old_rp:
            return workflow, resource_pool

    if 'func_SFS' not in resource_pool.keys():
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = calculate_sfs_workflow(workflow, resource_pool, config, name)
        if resource_pool == old_rp:
            return workflow, resource_pool

    fd = nipype_pipe_engine.Node(nipype_utility.Function(
        input_names=['in_file'], output_names=['out_file'],
        function=temporal_qc.fd_jenkinson), name='generate_FD_file{0}'.format(name))

    if 'mcflirt_rel_rms' in resource_pool.keys():
        fd.inputs.in_file = resource_pool['mcflirt_rel_rms']
    else:
        if len(resource_pool['func_coordinate_transformation']) == 2:
            node, out_file = resource_pool['func_coordinate_transformation']
            workflow.connect(node, out_file, fd, 'in_file')
        else:
            fd.inputs.in_file = \
                resource_pool['func_coordinate_transformation']

    spatial_epi = nipype_pipe_engine.Node(nipype_utility.Function(
        input_names=['mean_epi', 'func_brain_mask', 'direction', 'subject_id',
                     'session_id', 'scan_id', 'run_name', 'site_name',
                     'starter'],
        output_names=['qap'], function=qap_workflows_utils.qap_functional_spatial),
        name='qap_functional_spatial{0}'.format(name))

    # Subject infos
    if 'ghost_direction' not in config.keys():
        config['ghost_direction'] = 'y'

    spatial_epi.inputs.direction = config['ghost_direction']
    spatial_epi.inputs.subject_id = config['subject_id']
    spatial_epi.inputs.session_id = config['session_id']
    spatial_epi.inputs.scan_id = config['scan_id']
    spatial_epi.inputs.run_name = config['run_name']

    if 'site_name' in config.keys():
        spatial_epi.inputs.site_name = config['site_name']

    if len(resource_pool['func_mean']) == 2:
        node, out_file = resource_pool['func_mean']
        workflow.connect(node, out_file, spatial_epi, 'mean_epi')
    else:
        spatial_epi.inputs.mean_epi = resource_pool['func_mean']

    if len(resource_pool['func_brain_mask']) == 2:
        node, out_file = resource_pool['func_brain_mask']
        workflow.connect(node, out_file, spatial_epi, 'func_brain_mask')
    else:
        spatial_epi.inputs.func_brain_mask = \
            resource_pool['func_brain_mask']

    if config.get('write_report', False):
        plot = nipype_pipe_engine.Node(qap_viz.PlotMosaic(), name='func_plot_mosaic{0}'.format(name))
        plot.inputs.subject = config['subject_id']

        metadata = [config['session_id'], config['scan_id']]
        if 'site_name' in config.keys():
            metadata.append(config['site_name'])

        plot.inputs.metadata = metadata
        plot.inputs.title = 'Mean EPI'

        if len(resource_pool['func_mean']) == 2:
            node, out_file = resource_pool['func_mean']
            workflow.connect(node, out_file, plot, 'in_file')
        else:
            plot.inputs.in_file = resource_pool['func_mean']

        resource_pool['qap_mosaic'] = (plot, 'out_file')

    temporal = nipype_pipe_engine.Node(nipype_utility.Function(
        input_names=['func_timeseries', 'func_mean', 'func_brain_mask',
                     'bg_func_brain_mask', 'fd_file', 'motion_file', 'sfs',
                     'quality', 'outliers', 'oob_outliers', 'subject_id',
                     'session_id', 'scan_id', 'run_name', 'site_name',
                     'session_output_dir', 'starter'],
        output_names=['qap'],
        function=qap_workflows_utils.qap_functional_temporal),
        name='qap_functional_temporal{0}'.format(name))
    temporal.inputs.subject_id = config['subject_id']
    temporal.inputs.session_id = config['session_id']
    temporal.inputs.scan_id = config['scan_id']
    temporal.inputs.run_name = config['run_name']
    temporal.inputs.session_output_dir = config['output_directory']
    workflow.connect(fd, 'out_file', temporal, 'fd_file')

    if 'site_name' in config.keys():
        temporal.inputs.site_name = config['site_name']

    tqual = nipype_pipe_engine.Node(nipype_afni_preprocess.QualityIndex(), name='functional_quality_index%s' % name)
    workflow.connect(tqual, 'out_file', temporal, 'quality')

    gs_ts = nipype_pipe_engine.Node(nipype_utility.Function(input_names=["functional_file"],
                                                            output_names=["output"],
                                                            function=qap_workflows_utils.global_signal_time_series),
                                    name="global_signal_time_series{0}".format(name))
    
    outliers = nipype_pipe_engine.Node(nipype_afni_preprocess.OutlierCount(), name='func_outliers%s' % name)
    outliers.inputs.fraction = True
    outliers.inputs.save_outliers = False
    outliers.inputs.out_file = 'func_outliers%s' % name
    workflow.connect(outliers, 'out_file', temporal, 'outliers')

    oob_outliers = nipype_pipe_engine.Node(nipype_afni_preprocess.OutlierCount(), name='func_oob_outliers%s' % name)
    oob_outliers.inputs.fraction = True
    oob_outliers.inputs.save_outliers = False
    oob_outliers.inputs.out_file = 'func_oob_outliers%s' % name
    workflow.connect(oob_outliers, 'out_file', temporal, 'oob_outliers')

    if len(resource_pool['func_coordinate_transformation']) == 2:
        node, out_file = resource_pool['func_coordinate_transformation']
        workflow.connect(node, out_file, temporal, 'motion_file')
    else:
        temporal.inputs.motion_file = \
            resource_pool['func_coordinate_transformation']

    # func reorient (timeseries) -> QAP func temp
    if len(resource_pool['func_reorient']) == 2:
        node, out_file = resource_pool['func_reorient']
        workflow.connect(node, out_file, temporal, 'func_timeseries')
        workflow.connect(node, out_file, gs_ts, 'functional_file')
        workflow.connect(node, out_file, tqual, 'in_file')
        workflow.connect(node, out_file, outliers, 'in_file')
        workflow.connect(node, out_file, oob_outliers, 'in_file')
    else:
        qap_utils.check_input_resources(resource_pool, 'func_reorient')
        input_file = resource_pool['func_reorient']
        temporal.inputs.func_timeseries = input_file
        gs_ts.inputs.functional_file = input_file
        tqual.inputs.in_file = input_file
        outliers.inputs.in_file = input_file
        oob_outliers.inputs.in_file = input_file

    # func mean (one volume) -> QAP func temp
    if len(resource_pool['func_mean']) == 2:
        node, out_file = resource_pool['func_mean']
        workflow.connect(node, out_file, temporal, 'func_mean')
    else:
        qap_utils.check_input_resources(resource_pool, 'func_mean')
        input_file = resource_pool['func_mean']
        temporal.inputs.func_mean = input_file

    # functional brain mask -> QAP func temp
    if len(resource_pool['func_brain_mask']) == 2:
        node, out_file = resource_pool['func_brain_mask']
        workflow.connect(node, out_file, temporal, 'func_brain_mask')
        workflow.connect(node, out_file, outliers, 'mask')
    else:
        temporal.inputs.func_brain_mask = \
        outliers.inputs.mask = \
            resource_pool['func_brain_mask']

    # inverted functional brain mask -> QAP func temp
    if len(resource_pool['func_inverted_brain_mask']) == 2:
        node, out_file = resource_pool['func_inverted_brain_mask']
        workflow.connect(node, out_file, temporal, 'bg_func_brain_mask')
        workflow.connect(node, out_file, oob_outliers, 'mask')
    else:
        temporal.inputs.bg_func_brain_mask = \
        oob_outliers.inputs.mask = \
            resource_pool['func_inverted_brain_mask']

    # temporal STD -> QAP func temp
    if isinstance(resource_pool['func_SFS'], tuple):
        node, out_file = resource_pool['func_SFS']
        workflow.connect(node, out_file, temporal, 'sfs')
    else:
        temporal.inputs.sfs = resource_pool['func_SFS']

    # Write mosaic and FD plot

    out_dir = os.path.join(config['output_directory'], config["run_name"],
                           config["site_name"], config["subject_id"],
                           config["session_id"])
    out_json = os.path.join(out_dir, "{0}_{1}_{2}_qap-func.json".format(config["subject_id"], config["session_id"],
                                                                        config["scan_id"]))

    spatial_epi_to_json = nipype_pipe_engine.Node(nipype_utility.Function(
                                  input_names=["output_dict",
                                               "json_file"],
                                  output_names=["json_file"],
                                  function=qap_utils.write_json),
                               name="qap_functional_spatial_to_json{0}".format(name))
    spatial_epi_to_json.inputs.json_file = out_json

    temporal_to_json = nipype_pipe_engine.Node(nipype_utility.Function(
                                  input_names=["output_dict",
                                               "json_file"],
                                  output_names=["json_file"],
                                  function=qap_utils.write_json),
                              name="qap_functional_temporal_to_json{0}".format(name))
    temporal_to_json.inputs.json_file = out_json

    workflow.connect(spatial_epi, 'qap', spatial_epi_to_json, 'output_dict')
    workflow.connect(temporal, 'qap', temporal_to_json, 'output_dict')
    resource_pool['qap_functional'] = out_json

    if config['write_report']:

        metadata = [config['session_id'], config['scan_id']]
        if 'site_name' in config.keys():
            metadata += [config['site_name']]

        # todo: fix code to new qap

        out_ts_measures = os.path.join(out_dir, "{0}_{1}_{2}_timeseries-measures.png".format(config["subject_id"],
                                                                                             config["session_id"],
                                                                                             config["scan_id"]))

        out_cluster = os.path.join(out_dir, "{0}_{1}_{2}_grayplot-cluster.nii.gz".format(config["subject_id"],
                                                                                         config["session_id"],
                                                                                         config["scan_id"]))

        grayplot = nipype_pipe_engine.Node(qap_viz.GrayPlot(), name='grayplot{0}'.format(name))
        grayplot.inputs.subject = config['subject_id']
        grayplot.inputs.out_file = out_ts_measures
        grayplot.inputs.out_cluster = out_cluster
        id_string = "{0} {1} {2}".format(config["subject_id"], config["session_id"], config["scan_id"])
        grayplot.inputs.metadata = [id_string]
        workflow.connect(fd, 'out_file', grayplot, 'meanfd_file')
        dict_id = "{0} {1} {2}".format(config["subject_id"], config["session_id"], config["scan_id"])
        workflow.connect(temporal, ('qap', pick_dvars, dict_id), grayplot, 'dvars')    
        workflow.connect(gs_ts, 'output', grayplot, 'global_signal')
        resource_pool['timeseries_measures'] = (grayplot, 'out_file')
        resource_pool['grayplot-cluster'] = (grayplot, 'out_cluster')

        if len(resource_pool['func_reorient']) == 2:
            node, out_file = resource_pool['func_reorient']
            workflow.connect(node, out_file, grayplot, 'func_file')
        else:
            input_file = resource_pool['func_reorient']
            grayplot.inputs.func_file = input_file

        if len(resource_pool['func_brain_mask']) == 2:
            node, out_file = resource_pool['func_brain_mask']
            workflow.connect(node, out_file, grayplot, 'mask_file')

        else:   
            grayplot.inputs.mask_file = resource_pool['func_brain_mask']

    return workflow, resource_pool


def pick_dvars(qa, dvars_dict_id):
    print(qa[dvars_dict_id])
    return qa[dvars_dict_id]['metrics']['Standardized DVARS']