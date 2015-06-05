
def qap_mask_workflow(workflow, resource_pool, config):

    import os
    import sys

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    import nipype.interfaces.utility as util
    import nipype.interfaces.fsl.maths as fsl

    from workflows.qap_mask_utils import select_thresh, \
                                         slice_head_mask

    from workflows.workflow_utils import check_input_resources, \
                                         check_config_settings

  
    #check_input_resources(resource_pool, "anatomical_reorient")
    #check_input_resources(resource_pool, "ants_affine_xfm")
    check_config_settings(config, "template_skull_for_anat")


    if "ants_affine_xfm" not in resource_pool.keys():
        
        from workflows.anatomical_preproc import \
                                          ants_anatomical_linear_registration

        workflow, resource_pool = \
            ants_anatomical_linear_registration(workflow, resource_pool, 
                                                config)


    if "anatomical_reorient" not in resource_pool.keys():

        from workflows.anatomical_preproc import anatomical_reorient_workflow

        workflow, resource_pool = \
            anatomical_reorient_workflow(workflow, resource_pool, config)


    select_thresh = pe.Node(util.Function(input_names=['input_skull'],
                                          output_names=['thresh_out'],
                                          function=select_thresh),
                                          name='qap_headmask_select_thresh',
                                          iterfield=['input_skull'])


    mask_skull = pe.Node(interface=fsl.Threshold(), \
                         name='qap_headmask_thresh')
                         
    mask_skull.inputs.args = "-bin"



    dilate_node = pe.Node(interface=fsl.MathsCommand(), \
                          name="qap_headmask_dilate")

    dilate_node.inputs.args = "-dilM -dilM -dilM -dilM -dilM -dilM"


    erode_node = pe.Node(interface=fsl.MathsCommand(), \
                         name="qap_headmask_erode")

    erode_node.inputs.args = "-eroF -eroF -eroF -eroF -eroF -eroF"


    slice_head_mask = pe.Node(util.Function(input_names=['infile', 'transform', 'standard'],
                                            output_names=['outfile_path'],
                                            function=slice_head_mask),
                                          name="qap_headmask_slice_head_mask")


    combine_masks = pe.Node(interface=fsl.BinaryMaths(),
                            name="qap_headmask_combine_masks")

    combine_masks.inputs.operation = "add"
    combine_masks.inputs.args = "-bin"



    if len(resource_pool["anatomical_reorient"]) == 2:
        node, out_file = resource_pool["anatomical_reorient"]
        workflow.connect(node, out_file, select_thresh, 'input_skull')
        workflow.connect(node, out_file, mask_skull, 'in_file')
        workflow.connect(node, out_file, slice_head_mask, 'infile')
    else:
        select_thresh.inputs.input_skull = \
            resource_pool["anatomical_reorient"]
        mask_skull.inputs.in_file = \
            resource_pool["anatomical_reorient"]
        slice_head_mask.inputs.infile = \
            resource_pool["anatomical_reorient"]

     
    if len(resource_pool["ants_affine_xfm"]) == 2:
        node, out_file = resource_pool["ants_affine_xfm"]
        workflow.connect(node, out_file, slice_head_mask, 'transform')
    else:
        slice_head_mask.inputs.transform = resource_pool["ants_affine_xfm"]
                
    
    slice_head_mask.inputs.standard = config["template_skull_for_anat"]



    workflow.connect(select_thresh, 'thresh_out', mask_skull, 'thresh')

    workflow.connect(mask_skull, 'out_file', dilate_node, 'in_file')

    workflow.connect(dilate_node, 'out_file', erode_node, 'in_file')

    workflow.connect(erode_node, 'out_file', combine_masks, 'in_file')

    workflow.connect(slice_head_mask, 'outfile_path', \
                         combine_masks, 'operand_file')


    resource_pool["qap_head_mask"] = (combine_masks, 'out_file')


    return workflow, resource_pool



def qap_spatial_workflow(workflow, resource_pool, config):

    # resource pool should have:
    #     anatomical_reorient
    #     qap_head_mask
    #     anatomical_gm_mask
    #     anatomical_wm_mask
    #     anatomical_csf_mask

    import os
    import sys

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    import nipype.interfaces.utility as util
    
    from workflows.qap_measures_utils import qap_spatial, \
                                             append_to_csv


    if "qap_head_mask" not in resource_pool.keys():
    
        from workflows.qap_workflows import qap_mask_workflow

        workflow, resource_pool = \
            qap_mask_workflow(workflow, resource_pool, config)


    if ("anatomical_gm_mask" not in resource_pool.keys()) or \
           ("anatomical_wm_mask" not in resource_pool.keys()) or \
               ("anatomical_csf_mask" not in resource_pool.keys()):

        from workflows.anatomical_preproc import segmentation_workflow

        workflow, resource_pool = \
            segmentation_workflow(workflow, resource_pool, config)

    
    if "anatomical_reorient" not in resource_pool.keys():

        from workflows.anatomical_preproc import anatomical_reorient_workflow

        workflow, resource_pool = \
            anatomical_reorient_workflow(workflow, resource_pool, config)
    

    spatial = pe.Node(util.Function(input_names=['anatomical_reorient',
                                                 'head_mask_path',
                                                 'anatomical_gm_mask',
                                                 'anatomical_wm_mask',
                                                 'anatomical_csf_mask',
                                                 'subject_id',
                                                 'session_id',
                                                 'scan_id'],
                                    output_names=['qc'],
                                    function=qap_spatial),
                                    name='qap_spatial')
                                      
                                        
    spatial_to_csv = pe.Node(util.Function(input_names=['sub_qap_dict',
                                                        'outfile',
                                                        'append'],
                                           output_names=['outfile'],
                                           function=append_to_csv),
                                           name='qap_spatial_to_csv')


    spatial_to_main_csv = pe.Node(util.Function(input_names=['sub_qap_dict',
                                                             'outfile',
                                                             'append'],
                                                output_names=['outfile'],
                                                function=append_to_csv),
                                               name='qap_spatial_to_main_csv')
                                                                                

    if len(resource_pool["anatomical_reorient"]) == 2:
        node, out_file = resource_pool["anatomical_reorient"]
        workflow.connect(node, out_file, spatial, 'anatomical_reorient')
    else:
        spatial.inputs.anatomical_reorient = \
            resource_pool["anatomical_reorient"]


    if len(resource_pool["qap_head_mask"]) == 2:
        node, out_file = resource_pool["qap_head_mask"]
        workflow.connect(node, out_file, spatial, 'head_mask_path')
    else:
        spatial.inputs.head_mask_path = resource_pool["qap_head_mask"]

    
    if len(resource_pool["anatomical_gm_mask"]) == 2:
        node, out_file = resource_pool["anatomical_gm_mask"]
        workflow.connect(node, out_file, spatial, 'anatomical_gm_mask')
    else:
        spatial.inputs.anatomical_gm_mask = \
            resource_pool["anatomical_gm_mask"]


    if len(resource_pool["anatomical_wm_mask"]) == 2:
        node, out_file = resource_pool["anatomical_wm_mask"]
        workflow.connect(node, out_file, spatial, 'anatomical_wm_mask')
    else:
        spatial.inputs.anatomical_wm_mask = \
            resource_pool["anatomical_wm_mask"]
            
            
    if len(resource_pool["anatomical_csf_mask"]) == 2:
        node, out_file = resource_pool["anatomical_csf_mask"]
        workflow.connect(node, out_file, spatial, 'anatomical_csf_mask')
    else:
        spatial.inputs.anatomical_csf_mask = \
            resource_pool["anatomical_csf_mask"]
              
    
    # Subject infos 
     
    spatial.inputs.subject_id = config["subject_id"]
    spatial.inputs.session_id = config["session_id"]
    spatial.inputs.scan_id = config["scan_id"]
            
    workflow.connect(spatial, 'qc', spatial_to_csv, 'sub_qap_dict')

    spatial_to_csv.inputs.outfile = os.path.join(os.getcwd(), \
                                     "qap_spatial.csv")

    spatial_to_csv.inputs.append = False

            
    resource_pool["qap_spatial"] = (spatial_to_csv, 'outfile')


    # append to the main CSV

    workflow.connect(spatial, 'qc', spatial_to_main_csv, 'sub_qap_dict')

    spatial_to_main_csv.inputs.outfile = \
            os.path.join(config["output_directory"], "qap_spatial.csv")

    spatial_to_main_csv.inputs.append = True
    
    
    return workflow, resource_pool



def qap_spatial_epi_workflow(workflow, resource_pool, config):

    # resource pool should have:
    #     mean_functional
    #     functional_brain_mask

    import os
    import sys

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    import nipype.interfaces.utility as util
    
    from workflows.qap_measures_utils import qap_spatial_epi, \
                                             append_to_csv

    from workflows.workflow_utils import check_input_resources
  

    if "mean_functional" not in resource_pool.keys():

        from workflows.functional_preproc import mean_functional_workflow

        workflow, resource_pool = \
            mean_functional_workflow(workflow, resource_pool, config)


    if "functional_brain_mask" not in resource_pool.keys():

        from workflows.functional_preproc import \
                                              functional_brain_mask_workflow

        workflow, resource_pool = \
            functional_brain_workflow(workflow, resource_pool, config)


    spatial_epi = pe.Node(util.Function(input_names=['mean_epi',
                                                     'func_brain_mask',
                                                     'subject_id',
                                                     'session_id',
                                                     'scan_id'],
                                        output_names=['qc'],
                                        function=qap_spatial_epi),
                                        name='qap_spatial_epi')
                                         
                                           
    spatial_epi_to_csv = pe.Node(util.Function(input_names=['sub_qap_dict',
                                                            'outfile',
                                                            'append'],
                                               output_names=['outfile'],
                                               function=append_to_csv),
                                               name='qap_spatial_epi_to_csv')


    spatial_epi_to_main_csv = pe.Node(util.Function(input_names= \
                                               ['sub_qap_dict',
                                                'outfile',
                                                'append'],
                                               output_names=['outfile'],
                                               function=append_to_csv),
                                           name='qap_spatial_epi_to_main_csv')                                       
                                              

    if len(resource_pool["mean_functional"]) == 2:
        node, out_file = resource_pool["mean_functional"]
        workflow.connect(node, out_file, spatial_epi, 'mean_epi')
    else:
        spatial_epi.inputs.mean_epi = resource_pool["mean_functional"]


    if len(resource_pool["functional_brain_mask"]) == 2:
        node, out_file = resource_pool["functional_brain_mask"]
        workflow.connect(node, out_file, spatial_epi, 'func_brain_mask')
    else:
        spatial_epi.inputs.func_brain_mask = \
            resource_pool["functional_brain_mask"]
          
    
    # Subject infos 
     
    spatial_epi.inputs.subject_id = config["subject_id"]
    spatial_epi.inputs.session_id = config["session_id"]      
    spatial_epi.inputs.scan_id = config["scan_id"]        
    
    workflow.connect(spatial_epi, 'qc', spatial_epi_to_csv, 'sub_qap_dict')

    spatial_epi_to_csv.inputs.outfile = os.path.join(os.getcwd(), \
                                     "qap_spatial_epi.csv")

    spatial_epi_to_csv.inputs.append = False

 
    resource_pool["qap_spatial_epi"] = (spatial_epi_to_csv, 'outfile')


    # append to the main CSV

    workflow.connect(spatial_epi, 'qc',
                         spatial_epi_to_main_csv, 'sub_qap_dict')

    spatial_epi_to_main_csv.inputs.outfile = \
            os.path.join(config["output_directory"], \
                "qap_spatial_epi_%s.csv" % config["scan_id"])

    spatial_epi_to_main_csv.inputs.append = True
    
    
    return workflow, resource_pool



def qap_temporal_workflow(workflow, resource_pool, config):

    # resource pool should have:
    #     functional_brain_mask
    #     func_motion_correct
    #     coordinate_transformation

    import os
    import sys

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    import nipype.interfaces.utility as util
    
    from workflows.qap_measures_utils import qap_temporal, \
                                             append_to_csv
   
    
    if "mean_functional" not in resource_pool.keys():

        from workflows.functional_preproc import mean_functional_workflow

        workflow, resource_pool = \
            mean_functional_workflow(workflow, resource_pool, config)


    if ("func_motion_correct" not in resource_pool.keys()) or \
           ("coordinate_transformation" not in resource_pool.keys()):

        from workflows.functional_preproc import func_motion_correct_workflow

        workflow, resource_pool = \
            func_motion_correct_workflow(workflow, resource_pool, config)


    temporal = pe.Node(util.Function(input_names=['func_motion_correct',
                                                  'func_brain_mask',
                                                  'coord_xfm_matrix',
                                                  'subject_id',
                                                  'session_id',
                                                  'scan_id'],
                                     output_names=['qc'],
                                     function=qap_temporal),
                                     name='qap_temporal')


    temporal_to_csv = pe.Node(util.Function(input_names=['sub_qap_dict',
                                                         'outfile',
                                                         'append'],
                                            output_names=['outfile'],
                                            function=append_to_csv),
                                            name='qap_temporal_to_csv')  


    temporal_to_main_csv = pe.Node(util.Function(input_names=['sub_qap_dict',
                                                              'outfile',
                                                              'append'],
                                                 output_names=['outfile'],
                                                 function=append_to_csv),
                                              name='qap_temporal_to_main_csv')                                  
                                                                               
    
    if len(resource_pool["func_motion_correct"]) == 2:
        node, out_file = resource_pool["func_motion_correct"]
        workflow.connect(node, out_file, temporal, 'func_motion_correct')
    else:
        temporal.inputs.func_motion_correct = \
            resource_pool["func_motion_correct"]


    if len(resource_pool["functional_brain_mask"]) == 2:
        node, out_file = resource_pool["functional_brain_mask"]
        workflow.connect(node, out_file, temporal, 'func_brain_mask')
    else:
        temporal.inputs.func_brain_mask = \
            resource_pool["functional_brain_mask"]
            
            
    if len(resource_pool["coordinate_transformation"]) == 2:
        node, out_file = resource_pool["coordinate_transformation"]
        workflow.connect(node, out_file, temporal, 'coord_xfm_matrix')
    else:
        temporal.inputs.coord_xfm_matrix = \
            resource_pool["coordinate_transformation"]
            
    
    # Subject infos 
       
    temporal.inputs.subject_id = config["subject_id"]
    temporal.inputs.session_id = config["session_id"]
    temporal.inputs.scan_id = config["scan_id"]
                
    workflow.connect(temporal, 'qc', temporal_to_csv, 'sub_qap_dict')

    temporal_to_csv.inputs.outfile = os.path.join(os.getcwd(), \
                                      "qap_temporal.csv")

    temporal_to_csv.inputs.append = False
    
    resource_pool["qap_temporal"] = (temporal_to_csv, 'outfile')


    # append to the main CSV

    workflow.connect(temporal, 'qc', temporal_to_main_csv, 'sub_qap_dict')

    temporal_to_main_csv.inputs.outfile = \
            os.path.join(config["output_directory"], \
                "qap_temporal_%s.csv" % config["scan_id"])
    
    temporal_to_main_csv.inputs.append = True
    
    
    return workflow, resource_pool



