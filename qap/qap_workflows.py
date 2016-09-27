#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import os.path as op


def qap_mask_workflow(workflow, resource_pool, config, name="_"):

    import os
    import sys
    import copy
    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu
    from nipype.interfaces.afni import preprocess

    from qap_workflows_utils import create_expr_string, slice_head_mask
    from workflow_utils import check_input_resources, check_config_settings

    if 'template_skull_for_anat' not in config:
        config['template_skull_for_anat'] = Info.standard_image(
            'MNI152_T1_2mm.nii.gz')

    check_config_settings(config, 'template_skull_for_anat')

    if 'allineate_linear_xfm' not in resource_pool.keys():

        from anatomical_preproc import afni_anatomical_linear_registration
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            afni_anatomical_linear_registration(workflow, resource_pool,
                                                    config, name)

        if resource_pool == old_rp:
            return workflow, resource_pool

    if 'anatomical_reorient' not in resource_pool.keys():

        from anatomical_preproc import anatomical_reorient_workflow
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            anatomical_reorient_workflow(workflow, resource_pool, config, name)

        if resource_pool == old_rp:
            return workflow, resource_pool

    # find the clipping level for thresholding the head mask
    clip_level = pe.Node(interface=preprocess.ClipLevel(),
                         name='qap_headmask_clip_level%s' % name)

    # create the expression string for Calc in mask_skull node
    create_expr_string = pe.Node(niu.Function(
        input_names=['clip_level_value'],
        output_names=['expr_string'], function=create_expr_string),
        name='qap_headmask_create_expr_string%s' % name)

    workflow.connect(clip_level, 'clip_val', 
                         create_expr_string, 'clip_level_value')

    # let's create a binary mask of the skull image with that threshold
    mask_skull = pe.Node(interface=preprocess.Calc(),
                         name='qap_headmask_mask_skull%s' % name)
    mask_skull.inputs.outputtype = "NIFTI_GZ"

    workflow.connect(create_expr_string, 'expr_string', mask_skull, 'expr')


    dilate_erode = pe.Node(interface=preprocess.MaskTool(),
                           name='qap_headmask_mask_tool%s' % name)

    dilate_erode.inputs.dilate_inputs = "6 -6"
    dilate_erode.inputs.outputtype = "NIFTI_GZ"

    workflow.connect(mask_skull, 'out_file', dilate_erode, 'in_file')

    slice_head_mask = pe.Node(niu.Function(
        input_names=['infile', 'transform', 'standard'],
        output_names=['outfile_path'], function=slice_head_mask),
        name='qap_headmask_slice_head_mask%s' % name)

    combine_masks = pe.Node(interface=preprocess.Calc(), 
                            name='qap_headmask_combine_masks%s' % name)

    combine_masks.inputs.expr = "(a+b)-(a*b)"
    combine_masks.inputs.outputtype = "NIFTI_GZ"

    # subtract the slice mask from the original head mask to create a
    # skull-only mask for FG calculations
    subtract_mask = pe.Node(interface=preprocess.Calc(),
                            name='qap_headmask_subtract_masks%s' % name)

    subtract_mask.inputs.expr = "a-b"
    subtract_mask.inputs.outputtype = "NIFTI_GZ"

    if len(resource_pool['anatomical_reorient']) == 2:
        node, out_file = resource_pool['anatomical_reorient']
        workflow.connect([
            (node, clip_level,      [(out_file, 'in_file')]),
            (node, mask_skull,      [(out_file, 'in_file_a')]),
            (node, slice_head_mask, [(out_file, 'infile')])
        ])
    else:
        clip_level.inputs.in_file = resource_pool['anatomical_reorient']
        mask_skull.inputs.in_file_a = resource_pool['anatomical_reorient']
        slice_head_mask.inputs.infile = resource_pool['anatomical_reorient']

    if len(resource_pool['allineate_linear_xfm']) == 2:
        node, out_file = resource_pool['allineate_linear_xfm']
        workflow.connect(node, out_file, slice_head_mask, 'transform')
    else:
        slice_head_mask.inputs.transform = \
            resource_pool['allineate_linear_xfm']

    slice_head_mask.inputs.standard = config['template_skull_for_anat']

    workflow.connect([
        (dilate_erode, combine_masks, [('out_file', 'in_file_a')]),
        (slice_head_mask, combine_masks, [('outfile_path', 'in_file_b')])
    ])

    workflow.connect(dilate_erode, 'out_file', subtract_mask, 'in_file_a')
    workflow.connect(slice_head_mask, 'outfile_path', \
                         subtract_mask, 'in_file_b')

    resource_pool['qap_head_mask'] = (combine_masks, 'out_file')
    resource_pool['whole_head_mask'] = (dilate_erode, 'out_file')
    resource_pool['skull_only_mask'] = (subtract_mask, 'out_file')

    return workflow, resource_pool


def run_qap_mask(anatomical_reorient, allineate_out_xfm, template_skull,
                 out_dir=None, run=True):

    # stand-alone runner for anatomical reorient workflow

    import os
    import sys

    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    output = 'qap_head_mask'

    workflow = pe.Workflow(name='%s_workflow' % output)

    if not out_dir:
        out_dir = os.getcwd()

    workflow_dir = os.path.join(out_dir, "workflow_output", output)
    workflow.base_dir = workflow_dir

    resource_pool = {}
    config = {}
    num_cores_per_subject = 1

    resource_pool['anatomical_reorient'] = anatomical_reorient
    resource_pool['allineate_linear_xfm'] = allineate_out_xfm
    config['template_skull_for_anat'] = template_skull

    workflow, resource_pool = \
        qap_mask_workflow(workflow, resource_pool, config)

    ds = pe.Node(nio.DataSink(), name='datasink_%s' % output)
    ds.inputs.base_directory = workflow_dir

    node, out_file = resource_pool[output]
    workflow.connect(node, out_file, ds, output)

    node, out_file = resource_pool["skull_only_mask"]
    workflow.connect(node, out_file, ds, "skull_only_mask")

    if run:
        workflow.run(
            plugin='MultiProc', plugin_args={'n_procs': num_cores_per_subject})
        outpath = glob.glob(os.path.join(workflow_dir, output, '*'))[0]
        return outpath

    else:
        return workflow, workflow.base_dir


def add_header_to_qap_dict(in_file, qap_dict=None):

    import nibabel
    from qap.workflow_utils import raise_smart_exception

    try:
        img = nibabel.load(in_file)
        img_header = img.header
    except:
        err = "You may not have an up-to-date installation of the Python " \
              "Nibabel package.\nYour Nibabel version: %s" % \
              str(nb.__version__)
        raise_smart_exception(locals(),err)

    if not qap_dict:
        qap_dict = {}

    info_labels = ["descrip", "db_name", "bitpix", "slice_start", \
                   "scl_slope", "scl_inter", "slice_end", "slice_duration", \
                   "toffset", "quatern_b", "quatern_c", "quatern_d", \
                   "qoffset_x", "qoffset_y", "qoffset_z", "srow_x", "srow_y",\
                   "srow_z", "aux_file", "intent_name", "slice_code", \
                   "data_type", "qform_code", "sform_code"]


    for info_label in info_labels:
        try:
            qap_dict[info_label] = str(img_header[info_label])
        except:
            print "\n\n%s field not in NIFTI header of %s\n\n" % \
                  (info_label, in_file)
            pass

    try:
        pixdim = img_header['pixdim']
        qap_dict["pix_dimx"] = str(pixdim[1])
        qap_dict["pix_dimy"] = str(pixdim[2])
        qap_dict["pix_dimz"] = str(pixdim[3])
        qap_dict["tr"] = str(pixdim[4])
    except:
        print "\n\npix_dim/TR fields not in NIFTI header of %s\n\n" % in_file
        pass

    try:
        qap_dict["extensions"] = len(img.header.extensions.get_codes())
    except:
        print "\n\nExtensions not in NIFTI header of %s\n\n" % in_file
        pass

    return qap_dict


def qap_anatomical_spatial_workflow(workflow, resource_pool, config, name="_",
                                    report=False):

    # resource pool should have:
    #     anatomical_reorient
    #     qap_head_mask
    #     anatomical_gm_mask
    #     anatomical_wm_mask
    #     anatomical_csf_mask

    import os
    import sys
    import copy
    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu
    import nipype.algorithms.misc as nam
    from qap_workflows_utils import qap_anatomical_spatial
    from qap.viz.interfaces import PlotMosaic
    from workflow_utils import check_config_settings

    check_config_settings(config, "template_skull_for_anat")

    if 'qap_head_mask' not in resource_pool.keys():

        from qap_workflows import qap_mask_workflow
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            qap_mask_workflow(workflow, resource_pool, config, name)

        if resource_pool == old_rp:
            return workflow, resource_pool

    if ('anatomical_gm_mask' not in resource_pool.keys()) or \
            ('anatomical_wm_mask' not in resource_pool.keys()) or \
            ('anatomical_csf_mask' not in resource_pool.keys()):

        from anatomical_preproc import afni_segmentation_workflow
        old_rp = copy.copy(resource_pool)
        workflow, new_resource_pool = \
            afni_segmentation_workflow(workflow, resource_pool, config, name)

        if resource_pool == old_rp:
            return workflow, resource_pool

    if 'anatomical_reorient' not in resource_pool.keys():
        
        from anatomical_preproc import anatomical_reorient_workflow
        old_rp = copy.copy(resource_pool)
        workflow, new_resource_pool = \
            anatomical_reorient_workflow(workflow, resource_pool, config, name)

        if resource_pool == old_rp:
            return workflow, resource_pool

    spatial = pe.Node(niu.Function(
        input_names=['anatomical_reorient', 'qap_head_mask_path',
                     'whole_head_mask_path', 'skull_mask_path',
                     'anatomical_gm_mask', 'anatomical_wm_mask',
                     'anatomical_csf_mask', 'subject_id',
                     'session_id', 'scan_id', 'site_name',
                     'starter'],
        output_names=['qc'], function=qap_anatomical_spatial),
        name='qap_anatomical_spatial%s' % name)

    # Subject infos
    spatial.inputs.subject_id = config['subject_id']
    spatial.inputs.session_id = config['session_id']
    spatial.inputs.scan_id = config['scan_id']

    node, out_file = resource_pool['starter']
    workflow.connect(node, out_file, spatial, 'starter')

    if 'site_name' in resource_pool.keys():
        spatial.inputs.site_name = resource_pool['site_name']

    if len(resource_pool['anatomical_reorient']) == 2:
        node, out_file = resource_pool['anatomical_reorient']
        workflow.connect(node, out_file, spatial, 'anatomical_reorient')
    else:
        spatial.inputs.anatomical_reorient = \
            resource_pool['anatomical_reorient']

    if len(resource_pool['qap_head_mask']) == 2:
        node, out_file = resource_pool['qap_head_mask']
        workflow.connect(node, out_file, spatial, 'qap_head_mask_path')
    else:
        spatial.inputs.qap_head_mask_path = resource_pool['qap_head_mask']

    if len(resource_pool['whole_head_mask']) == 2:
        node, out_file = resource_pool['whole_head_mask']
        workflow.connect(node, out_file, spatial, 'whole_head_mask_path')
    else:
        spatial.inputs.whole_head_mask_path = resource_pool['whole_head_mask']

    if len(resource_pool['skull_only_mask']) == 2:
        node, out_file = resource_pool['skull_only_mask']
        workflow.connect(node, out_file, spatial, 'skull_mask_path')
    else:
        spatial.inputs.skull_mask_path = resource_pool['skull_only_mask']

    if len(resource_pool['anatomical_gm_mask']) == 2:
        node, out_file = resource_pool['anatomical_gm_mask']
        workflow.connect(node, out_file, spatial, 'anatomical_gm_mask')
    else:
        spatial.inputs.anatomical_gm_mask = \
            resource_pool['anatomical_gm_mask']

    if len(resource_pool['anatomical_wm_mask']) == 2:
        node, out_file = resource_pool['anatomical_wm_mask']
        workflow.connect(node, out_file, spatial, 'anatomical_wm_mask')
    else:
        spatial.inputs.anatomical_wm_mask = \
            resource_pool['anatomical_wm_mask']

    if len(resource_pool['anatomical_csf_mask']) == 2:
        node, out_file = resource_pool['anatomical_csf_mask']
        workflow.connect(node, out_file, spatial, 'anatomical_csf_mask')
    else:
        spatial.inputs.anatomical_csf_mask = \
            resource_pool['anatomical_csf_mask']

    if config.get('write_report', False):
        plot = pe.Node(PlotMosaic(), name='plot_mosaic%s' % name)
        plot.inputs.subject = config['subject_id']

        metadata = [config['session_id'], config['scan_id']]
        if 'site_name' in config.keys():
            metadata.append(config['site_name'])

        plot.inputs.metadata = metadata
        plot.inputs.title = 'Anatomical reoriented'

        if len(resource_pool['anatomical_reorient']) == 2:
            node, out_file = resource_pool['anatomical_reorient']
            workflow.connect(node, out_file, plot, 'in_file')
        else:
            plot.inputs.in_file = resource_pool['anatomical_reorient']

        # Enable this if we want masks
        # if len(resource_pool['qap_head_mask']) == 2:
        #     node, out_file = resource_pool['qap_head_mask']
        #     workflow.connect(node, out_file, plot, 'in_mask')
        # else:
        #     plot.inputs.in_mask = resource_pool['qap_head_mask']
        
        resource_pool['qap_mosaic'] = (plot, 'out_file')

    add_header = pe.Node(niu.Function(
                             input_names=['in_file', 'qap_dict'],
                             output_names=['qap_dict'],
                             function=add_header_to_qap_dict),
                     name="add_header_to_anatomical_spatial_csv%s" % name)

    if len(resource_pool['anatomical_reorient']) == 2:
        node, out_file = resource_pool['anatomical_reorient']
        workflow.connect(node, out_file, add_header, 'in_file')
    else:
        add_header.inputs.in_file = resource_pool['anatomical_reorient']

    out_csv = op.join(config['output_directory'], 'qap_anatomical_spatial.csv')
    spatial_to_csv = pe.Node(
        nam.AddCSVRow(in_file=out_csv),
        name='qap_anatomical_spatial_to_csv%s' % name)

    workflow.connect(spatial, 'qc', add_header, 'qap_dict')
    workflow.connect(add_header, 'qap_dict', spatial_to_csv, '_outputs')
    resource_pool['qap_anatomical_spatial'] = (spatial_to_csv, 'csv_file')

    return workflow, resource_pool


def run_single_qap_anatomical_spatial(
        anatomical_reorient, qap_head_mask, anatomical_csf_mask,
        anatomical_gm_mask, anatomical_wm_mask, subject_id, session_id=None,
        scan_id=None, site_name=None, out_dir=None, run=True):

    # stand-alone runner for anatomical spatial QAP workflow

    import os
    import sys
    import glob
    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    output = 'qap_anatomical_spatial'
    workflow = pe.Workflow(name='%s_workflow' % output)

    if not out_dir:
        out_dir = os.getcwd()

    workflow_dir = os.path.join(out_dir, "workflow_output", output)
    workflow.base_dir = workflow_dir

    num_cores_per_subject = 1
    resource_pool = {
        'anatomical_reorient': anatomical_reorient,
        'qap_head_mask': qap_head_mask,
        'anatomical_csf_mask': anatomical_csf_mask,
        'anatomical_gm_mask': anatomical_gm_mask,
        'anatomical_wm_mask': anatomical_wm_mask
    }

    config = {
        'subject_id': subject_id,
        'session_id': session_id,
        'scan_id': scan_id
    }

    if site_name:
        config['site_name'] = site_name

    workflow, resource_pool = \
        qap_anatomical_spatial_workflow(workflow, resource_pool, config)

    ds = pe.Node(nio.DataSink(), name='datasink_%s' % output)
    ds.inputs.base_directory = workflow_dir

    node, out_file = resource_pool[output]

    workflow.connect(node, out_file, ds, output)

    if run:
        workflow.run(
            plugin='MultiProc', plugin_args={'n_procs': num_cores_per_subject})
        outpath = glob.glob(os.path.join(workflow_dir, output, '*'))[0]
        return outpath

    else:
        return workflow, workflow.base_dir


def run_whole_single_qap_anatomical_spatial(
        anatomical_scan, template_head, subject_id, session_id=None,
        scan_id=None, site_name=None, out_dir=None, run=True):

    # stand-alone runner for anatomical spatial QAP workflow

    import os
    import sys
    import glob
    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe

    from qap import cli

    output = 'qap_anatomical_spatial'
    workflow = pe.Workflow(name='%s_workflow' % output)

    if not out_dir:
        out_dir = os.getcwd()

    if site_name:
        workflow_dir = os.path.join(out_dir, "workflow_output", output, \
            site_name, subject_id)
    else:
        workflow_dir = os.path.join(out_dir, "workflow_output", output, \
            subject_id)

    if session_id:
        workflow_dir = os.path.join(workflow_dir, session_id)

    if scan_id:
        workflow_dir = os.path.join(workflow_dir, scan_id)

    workflow.base_dir = workflow_dir

    num_cores_per_subject = 1
    resource_pool = {
        'anatomical_scan': anatomical_scan
    }

    config = {
        'template_skull_for_anat': template_head,
        'subject_id': subject_id,
        'session_id': session_id,
        'scan_id': scan_id,
        'output_directory': workflow_dir
    }

    if site_name:
        config['site_name'] = site_name

    # create the one node all participants will start from
    starter_node = pe.Node(niu.Function(input_names=['starter'], 
                                        output_names=['starter'], 
                                        function=cli.starter_node_func),
                           name='starter_node')

    # set a dummy variable
    starter_node.inputs.starter = ""

    resource_pool["starter"] = (starter_node, 'starter')

    workflow, resource_pool = \
        qap_anatomical_spatial_workflow(workflow, resource_pool, config)

    ds = pe.Node(nio.DataSink(), name='datasink_%s' % output)
    ds.inputs.base_directory = workflow_dir

    node, out_file = resource_pool[output]

    workflow.connect(node, out_file, ds, output)

    if run:
        workflow.run(
            plugin='MultiProc', plugin_args={'n_procs': num_cores_per_subject})
        outpath = glob.glob(os.path.join(workflow_dir, output, '*'))[0]
        return outpath

    else:
        return workflow, workflow.base_dir


def qap_functional_spatial_workflow(workflow, resource_pool, config, name="_"):

    # resource pool should have:
    #     mean_functional
    #     functional_brain_mask

    import os
    import sys
    import copy
    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    import nipype.algorithms.misc as nam
    import nipype.interfaces.utility as niu
    import nipype.algorithms.misc as nam

    from qap_workflows_utils import qap_functional_spatial
    from qap.viz.interfaces import PlotMosaic

    from workflow_utils import check_input_resources

    if 'mean_functional' not in resource_pool.keys():
        from functional_preproc import mean_functional_workflow
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            mean_functional_workflow(workflow, resource_pool, config, name)
        if resource_pool == old_rp:
            return workflow, resource_pool

    if 'functional_brain_mask' not in resource_pool.keys():
        from functional_preproc import functional_brain_mask_workflow
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            functional_brain_mask_workflow(workflow, resource_pool, config, name)
        if resource_pool == old_rp:
            return workflow, resource_pool

    spatial_epi = pe.Node(niu.Function(
        input_names=['mean_epi', 'func_brain_mask', 'direction', 'subject_id',
                     'session_id', 'scan_id', 'site_name', 'starter'],
        output_names=['qc'], function=qap_functional_spatial),
        name='qap_functional_spatial%s' % name)

    # Subject infos
    if 'ghost_direction' not in config.keys():
        config['ghost_direction'] = 'y'

    spatial_epi.inputs.direction = config['ghost_direction']
    spatial_epi.inputs.subject_id = config['subject_id']
    spatial_epi.inputs.session_id = config['session_id']
    spatial_epi.inputs.scan_id = config['scan_id']

    if 'site_name' in config.keys():
        spatial_epi.inputs.site_name = config['site_name']

    if len(resource_pool['mean_functional']) == 2:
        node, out_file = resource_pool['mean_functional']
        workflow.connect(node, out_file, spatial_epi, 'mean_epi')
    else:
        spatial_epi.inputs.mean_epi = resource_pool['mean_functional']

    if len(resource_pool['functional_brain_mask']) == 2:
        node, out_file = resource_pool['functional_brain_mask']
        workflow.connect(node, out_file, spatial_epi, 'func_brain_mask')
    else:
        spatial_epi.inputs.func_brain_mask = \
            resource_pool['functional_brain_mask']

    if config.get('write_report', False):
        plot = pe.Node(PlotMosaic(), name='plot_mosaic%s' % name)
        plot.inputs.subject = config['subject_id']

        metadata = [config['session_id'], config['scan_id']]
        if 'site_name' in config.keys():
            metadata.append(config['site_name'])

        plot.inputs.metadata = metadata
        plot.inputs.title = 'Mean EPI'

        if len(resource_pool['mean_functional']) == 2:
            node, out_file = resource_pool['mean_functional']
            workflow.connect(node, out_file, plot, 'in_file')
        else:
            plot.inputs.in_file = resource_pool['mean_functional']

        # Enable this if we want masks
        # if len(resource_pool['functional_brain_mask']) == 2:
        #     node, out_file = resource_pool['functional_brain_mask']
        #     workflow.connect(node, out_file, plot, 'in_mask')
        # else:
        #     plot.inputs.in_mask = resource_pool['functional_brain_mask']
        resource_pool['qap_mosaic'] = (plot, 'out_file')

    add_header = pe.Node(niu.Function(
                             input_names=['in_file', 'qap_dict'],
                             output_names=['qap_dict'],
                             function=add_header_to_qap_dict),
                     name="add_header_to_functional_spatial_csv%s" % name)

    if len(resource_pool['mean_functional']) == 2:
        node, out_file = resource_pool['mean_functional']
        workflow.connect(node, out_file, add_header, 'in_file')
    else:
        add_header.inputs.in_file = resource_pool['mean_functional']

    out_csv = op.join(
        config['output_directory'], 'qap_functional_spatial.csv')
    spatial_epi_to_csv = pe.Node(
        nam.AddCSVRow(in_file=out_csv),
        name='qap_functional_spatial_to_csv%s' % name)

    workflow.connect(spatial_epi, 'qc', add_header, 'qap_dict')
    workflow.connect(add_header, 'qap_dict', spatial_epi_to_csv, '_outputs')
    
    resource_pool['qap_functional_spatial'] = (spatial_epi_to_csv, 'csv_file')

    return workflow, resource_pool


def run_single_qap_functional_spatial(
        mean_functional, functional_brain_mask, subject_id, session_id,
        scan_id, site_name=None, ghost_direction=None, out_dir=None,
        run=True):

    # stand-alone runner for functional spatial QAP workflow
    import os
    import sys
    import glob
    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    output = 'qap_functional_spatial'
    workflow = pe.Workflow(name='%s_workflow' % output)

    if not out_dir:
        out_dir = os.getcwd()

    workflow_dir = os.path.join(out_dir, "workflow_output", output)
    workflow.base_dir = workflow_dir

    resource_pool = {}
    config = {}
    num_cores_per_subject = 1

    resource_pool['mean_functional'] = mean_functional
    resource_pool['functional_brain_mask'] = functional_brain_mask

    config['subject_id'] = subject_id
    config['session_id'] = session_id
    config['scan_id'] = scan_id

    if site_name:
        config['site_name'] = site_name

    if ghost_direction:
        config['ghost_direction'] = ghost_direction

    workflow, resource_pool = \
        qap_functional_spatial_workflow(workflow, resource_pool, config)

    ds = pe.Node(nio.DataSink(), name='datasink_%s' % output)
    ds.inputs.base_directory = workflow_dir

    node, out_file = resource_pool[output]

    workflow.connect(node, out_file, ds, output)

    if run:
        workflow.run(
            plugin='MultiProc', plugin_args={'n_procs': num_cores_per_subject})
        outpath = glob.glob(os.path.join(workflow_dir, output, '*'))[0]
        return outpath

    else:
        return workflow, workflow.base_dir


def run_whole_single_qap_functional_spatial(
        functional_scan, subject_id, session_id=None, scan_id=None,
        site_name=None, out_dir=None, run=True):

    # stand-alone runner for functional spatial QAP workflow

    import os
    import sys
    import glob
    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe

    from qap import cli

    output = 'qap_functional_spatial'
    workflow = pe.Workflow(name='%s_workflow' % output)

    if not out_dir:
        out_dir = os.getcwd()

    if site_name != None:
        workflow_dir = os.path.join(out_dir, "workflow_output", output, \
            site_name, subject_id)
    else:
        workflow_dir = os.path.join(out_dir, "workflow_output", output, \
            subject_id)

    if session_id != None:
        workflow_dir = os.path.join(workflow_dir, session_id)

    if scan_id != None:
        workflow_dir = os.path.join(workflow_dir, scan_id)

    workflow.base_dir = workflow_dir

    num_cores_per_subject = 1
    resource_pool = {
        'functional_scan': functional_scan
    }

    config = {
        'subject_id': subject_id,
        'session_id': session_id,
        'scan_id': scan_id,
        'output_directory': workflow_dir
    }

    if site_name:
        config['site_name'] = site_name

    # create the one node all participants will start from
    starter_node = pe.Node(niu.Function(input_names=['starter'], 
                                        output_names=['starter'], 
                                        function=cli.starter_node_func),
                           name='starter_node')

    # set a dummy variable
    starter_node.inputs.starter = ""

    resource_pool["starter"] = (starter_node, 'starter')

    workflow, resource_pool = \
        qap_functional_spatial_workflow(workflow, resource_pool, config)

    ds = pe.Node(nio.DataSink(), name='datasink_%s' % output)
    ds.inputs.base_directory = workflow_dir

    node, out_file = resource_pool[output]

    workflow.connect(node, out_file, ds, output)

    if run:
        workflow.run(
            plugin='MultiProc', plugin_args={'n_procs': num_cores_per_subject})
        outpath = glob.glob(os.path.join(workflow_dir, output, '*'))[0]
        return outpath

    else:
        return workflow, workflow.base_dir


def qap_functional_temporal_workflow(workflow, resource_pool, config, name="_"):

    # resource pool should have:
    #     functional_brain_mask
    #     func_motion_correct
    #     coordinate_transformation

    import os
    import sys
    import copy
    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu
    import nipype.algorithms.misc as nam

    from qap_workflows_utils import qap_functional_temporal
    from temporal_qc import fd_jenkinson
    from qap.viz.interfaces import PlotMosaic, PlotFD

    def _getfirst(inlist):
        if isinstance(inlist, list):
            return inlist[0]
        return inlist

    # ensures functional_brain_mask is created as well
    if 'inverted_functional_brain_mask' not in resource_pool.keys():
        from functional_preproc import invert_functional_brain_mask_workflow
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            invert_functional_brain_mask_workflow(workflow, resource_pool, config, name)
        if resource_pool == old_rp:
            return workflow, resource_pool

    if ('func_motion_correct' not in resource_pool.keys()) or \
        ('coordinate_transformation' not in resource_pool.keys() and
            'mcflirt_rel_rms' not in resource_pool.keys()):
        from functional_preproc import func_motion_correct_workflow
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            func_motion_correct_workflow(workflow, resource_pool, config, name)
        if resource_pool == old_rp:
            return workflow, resource_pool

    fd = pe.Node(niu.Function(
        input_names=['in_file'], output_names=['out_file'],
        function=fd_jenkinson), name='generate_FD_file%s' % name)

    if 'mcflirt_rel_rms' in resource_pool.keys():
        fd.inputs.in_file = resource_pool['mcflirt_rel_rms']
    else:
        if len(resource_pool['coordinate_transformation']) == 2:
            node, out_file = resource_pool['coordinate_transformation']
            workflow.connect(node, out_file, fd, 'in_file')
        else:
            fd.inputs.in_file = resource_pool['coordinate_transformation']

    temporal = pe.Node(niu.Function(
        input_names=['func_timeseries', 'func_brain_mask',
                     'bg_func_brain_mask', 'fd_file', 'subject_id',
                     'session_id', 'scan_id', 'site_name', 'starter'],
        output_names=['qc'],
        function=qap_functional_temporal),
        name='qap_functional_temporal%s' % name)
    temporal.inputs.subject_id = config['subject_id']
    temporal.inputs.session_id = config['session_id']
    temporal.inputs.scan_id = config['scan_id']
    workflow.connect(fd, 'out_file', temporal, 'fd_file')

    if 'site_name' in config.keys():
        temporal.inputs.site_name = config['site_name']

    # func reorient (timeseries) -> QAP func temp
    if len(resource_pool['func_reorient']) == 2:
        node, out_file = resource_pool['func_reorient']
        workflow.connect(node, out_file, temporal, 'func_timeseries')
    else:
        from workflow_utils import check_input_resources
        check_input_resources(resource_pool, 'func_reorient')
        input_file = resource_pool['func_reorient']
        temporal.inputs.func_timeseries = input_file

    # functional brain mask -> QAP func temp
    if len(resource_pool['functional_brain_mask']) == 2:
        node, out_file = resource_pool['functional_brain_mask']
        workflow.connect(node, out_file, temporal, 'func_brain_mask')
    else:
        temporal.inputs.func_brain_mask = \
            resource_pool['functional_brain_mask']

    # inverted functional brain mask -> QAP func temp
    if len(resource_pool['inverted_functional_brain_mask']) == 2:
        node, out_file = resource_pool['inverted_functional_brain_mask']
        workflow.connect(node, out_file, temporal, 'bg_func_brain_mask')
    else:
        temporal.inputs.bg_func_brain_mask = \
            resource_pool['inverted_functional_brain_mask']

    # Write mosaic and FD plot
    if config.get('write_report', False):
        #plot = pe.Node(PlotMosaic(), name='plot_mosaic')
        #plot.inputs.subject = config['subject_id']

        metadata = [config['session_id'], config['scan_id']]
        if 'site_name' in config.keys():
            metadata.append(config['site_name'])

        #plot.inputs.metadata = metadata
        #plot.inputs.title = 'tSNR volume'
        #workflow.connect(tsnr, 'tsnr_file', plot, 'in_file')

        # Enable this if we want masks
        # if len(resource_pool['functional_brain_mask']) == 2:
        #     node, out_file = resource_pool['functional_brain_mask']
        #     workflow.connect(node, out_file, plot, 'in_mask')
        # else:
        #     plot.inputs.in_mask = resource_pool['functional_brain_mask']
        #resource_pool['qap_mosaic'] = (plot, 'out_file')

        fdplot = pe.Node(PlotFD(), name='plot_fd%s' % name)
        fdplot.inputs.subject = config['subject_id']
        fdplot.inputs.metadata = metadata
        workflow.connect(fd, 'out_file', fdplot, 'in_file')
        resource_pool['qap_fd'] = (fdplot, 'out_file')

    add_header = pe.Node(niu.Function(
                             input_names=['in_file', 'qap_dict'],
                             output_names=['qap_dict'],
                             function=add_header_to_qap_dict),
                             name="add_header_to_functional_temporal_csv%s" % name)

    if len(resource_pool['func_reorient']) == 2:
        node, out_file = resource_pool['func_reorient']
        workflow.connect(node, out_file, add_header, 'in_file')
    else:
        add_header.inputs.in_file = resource_pool['func_reorient']

    out_csv = op.join(
        config['output_directory'], 'qap_functional_temporal.csv')
    temporal_to_csv = pe.Node(
        nam.AddCSVRow(in_file=out_csv),
        name='qap_functional_temporal_to_csv%s' % name)

    workflow.connect(temporal, 'qc', add_header, 'qap_dict')
    workflow.connect(add_header, 'qap_dict', temporal_to_csv, '_outputs')

    resource_pool['qap_functional_temporal'] = (temporal_to_csv, 'csv_file')

    return workflow, resource_pool


def run_single_qap_functional_temporal(func_reorient, functional_brain_mask,
                                       subject_id, session_id, scan_id,
                                       site_name=None, mcflirt_rel_rms=None,
                                       coordinate_transformation=None,
                                       out_dir=None, run=True):

    # stand-alone runner for functional temporal QAP workflow

    import os
    import sys
    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    output = 'qap_functional_temporal'

    workflow = pe.Workflow(name='%s_workflow' % output)

    if not out_dir:
        out_dir = os.getcwd()

    workflow_dir = os.path.join(out_dir, "workflow_output", output)
    workflow.base_dir = workflow_dir

    resource_pool = {}
    config = {}
    num_cores_per_subject = 1

    resource_pool['func_reorient'] = func_reorient
    resource_pool['functional_brain_mask'] = functional_brain_mask

    if mcflirt_rel_rms:
        resource_pool['mcflirt_rel_rms'] = mcflirt_rel_rms
    elif coordinate_transformation:
        resource_pool['coordinate_transformation'] = coordinate_transformation

    config['subject_id'] = subject_id
    config['session_id'] = session_id
    config['scan_id'] = scan_id

    if site_name:
        config['site_name'] = site_name

    workflow, resource_pool = \
        qap_functional_temporal_workflow(workflow, resource_pool, config)

    ds = pe.Node(nio.DataSink(), name='datasink_%s' % output)
    ds.inputs.base_directory = workflow_dir

    node, out_file = resource_pool[output]

    workflow.connect(node, out_file, ds, output)

    if run:
        workflow.run(
            plugin='MultiProc', plugin_args={'n_procs': num_cores_per_subject})
        outpath = glob.glob(os.path.join(workflow_dir, output, '*'))[0]
        return outpath
    else:
        return workflow, workflow.base_dir



def run_whole_single_qap_functional_temporal(
        functional_scan, subject_id, session_id=None, scan_id=None,
        site_name=None, out_dir=None, run=True):

    # stand-alone runner for functional temporal QAP workflow

    import os
    import sys
    import glob
    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe

    from qap import cli

    output = 'qap_functional_temporal'
    workflow = pe.Workflow(name='%s_workflow' % output)

    if not out_dir:
        out_dir = os.getcwd()

    if site_name:
        workflow_dir = os.path.join(out_dir, "workflow_output", output, \
            site_name, subject_id)
    else:
        workflow_dir = os.path.join(out_dir, "workflow_output", output, \
            subject_id)

    if session_id:
        workflow_dir = os.path.join(workflow_dir, session_id)

    if scan_id:
        workflow_dir = os.path.join(workflow_dir, scan_id)

    workflow.base_dir = workflow_dir

    num_cores_per_subject = 1
    resource_pool = {
        'functional_scan': functional_scan
    }

    config = {
        'subject_id': subject_id,
        'session_id': session_id,
        'scan_id': scan_id,
        'output_directory': workflow_dir
    }

    if site_name:
        config['site_name'] = site_name

    # create the one node all participants will start from
    starter_node = pe.Node(niu.Function(input_names=['starter'], 
                                        output_names=['starter'], 
                                        function=cli.starter_node_func),
                           name='starter_node')

    # set a dummy variable
    starter_node.inputs.starter = ""

    resource_pool["starter"] = (starter_node, 'starter')

    workflow, resource_pool = \
        qap_functional_temporal_workflow(workflow, resource_pool, config)

    ds = pe.Node(nio.DataSink(), name='datasink_%s' % output)
    ds.inputs.base_directory = workflow_dir

    node, out_file = resource_pool[output]

    workflow.connect(node, out_file, ds, output)

    if run:
        workflow.run(
            plugin='MultiProc', plugin_args={'n_procs': num_cores_per_subject})
        outpath = glob.glob(os.path.join(workflow_dir, output, '*'))[0]
        return outpath

    else:
        return workflow, workflow.base_dir
