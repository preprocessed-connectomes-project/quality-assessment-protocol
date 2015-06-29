

test_sub_dir = "test_data/1019436/session_1"


def test_run_anatomical_reorient():

    import os
    import commands

    import pkg_resources as p

    from qap.anatomical_preproc import run_anatomical_reorient


    if "anatomical_reorient" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)


    anat_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_1", \
                                    "anatomical_scan", \
                                    "mprage.nii.gz"))

    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "anat_1", \
                                  "anatomical_reorient", \
                                  "mprage_resample.nii.gz"))

    # run the workflow
    output = run_anatomical_reorient(anat_scan)

    # calculate the correlation
    cmd = "3ddot -demean %s %s" % (output, ref_out)
    correlation = commands.getoutput(cmd)

    os.system("rm -R anatomical_reorient")


    assert correlation > 0.98
    
    
    
def test_run_anatomical_skullstrip():

    import os
    import commands
    
    import pkg_resources as p

    from qap.anatomical_preproc import run_anatomical_skullstrip


    if "anatomical_skullstrip" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)


    anat_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "anat_1", \
                                        "anatomical_reorient", \
                                        "mprage_resample.nii.gz"))

    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "anat_1", \
                                  "anatomical_brain", \
                                  "mprage_resample_calc.nii.gz"))

    # run the workflow
    output = run_anatomical_skullstrip(anat_reorient)

    # calculate the correlation
    cmd = "3ddot -demean %s %s" % (output, ref_out)
    correlation = commands.getoutput(cmd)

    os.system("rm -R anatomical_skullstrip")


    assert correlation > 0.98



def test_run_ants_anatomical_linear_registration():

    import os
    import commands
    
    import pkg_resources as p

    from qap.anatomical_preproc import run_ants_anatomical_linear_registration


    if "ants_anatomical_linear_registration" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)


    anat_brain = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_1", \
                                     "anatomical_brain", \
                                     "mprage_resample_calc.nii.gz"))

    template_brain = p.resource_filename("qap", os.path.join("test_data", \
                                         "MNI152_T1_2mm_brain.nii.gz"))

    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "anat_1", \
                                  "ants_linear_warped_image", \
                                  "transform_Warped.nii.gz"))

    # run the workflow
    output = run_ants_anatomical_linear_registration(anat_brain, \
                                                         template_brain)

    # calculate the correlation
    cmd = "3ddot -demean %s %s" % (output, ref_out)
    correlation = commands.getoutput(cmd)

    os.system("rm -R ants_anatomical_linear_registration")


    assert correlation > 0.98
    
    
    
def test_run_segmentation_workflow():

    import os
    import commands
    
    import pkg_resources as p

    from qap.anatomical_preproc import run_segmentation_workflow


    if "segmentation" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)


    anat_brain = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_1", \
                                     "anatomical_brain", \
                                     "mprage_resample_calc.nii.gz"))

    ref_csf_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                      "anat_1", \
                                      "anatomical_csf_mask", \
                                      "segment_seg_0.nii.gz"))
                                  
    ref_gm_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_1", \
                                     "anatomical_gm_mask", \
                                     "segment_seg_1.nii.gz"))
                                  
    ref_wm_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_1", \
                                     "anatomical_wm_mask", \
                                     "segment_seg_2.nii.gz"))

    # run the workflow
    output = run_segmentation_workflow(anat_brain, 0.98, 0.7, 0.98)

    ref_list = [ref_csf_out, ref_gm_out, ref_wm_out]

    correlation_count = 0

    # calculate the correlation
    for out, ref in zip(output, ref_list):
    
        cmd = "3ddot -demean %s %s" % (out, ref)
        correlation = commands.getoutput(cmd)
        
        if correlation > 0.98:
            correlation_count += 1

    os.system("rm -R segmentation")


    assert correlation_count == 3



def run_all_tests_anatomical_preproc():

    test_run_anatomical_reorient()
    test_run_anatomical_skullstrip()
    test_run_ants_anatomical_linear_registration()
    test_run_segmentation_workflow()



