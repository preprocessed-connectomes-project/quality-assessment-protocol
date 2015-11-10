

test_sub_dir = "test_data/1019436/session_1"


def test_run_func_motion_correct_no_slice_time():

    import os
    import nibabel as nb

    import pkg_resources as p

    from qap.functional_preproc import run_func_motion_correct


    if "func_motion_correct" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)

    
    func_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "functional_scan", \
                                    "rest.nii.gz"))
                                    
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "rest_1", \
                                  "func_motion_correct", \
                                  "rest_calc_tshift_resample_volreg.nii.gz"))
                                    
    # run the workflow
    output = run_func_motion_correct(func_scan, 0, "End", False)

    # make the correlation
    ref_out_data = nb.load(ref_out).get_data()  
    output_data = nb.load(output).get_data()
    
    os.system("rm -R func_motion_correct")

    # create a vector of True and False values
    bool_vector = ref_out_data == output_data

    assert bool_vector.all()



def test_run_func_motion_correct_slice_time():

    import os
    import nibabel as nb

    import pkg_resources as p

    from qap.functional_preproc import run_func_motion_correct


    if "func_motion_correct" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)

    
    func_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "functional_scan", \
                                    "rest.nii.gz"))
    
    ''' NEED A SLICE TIME CORRECTED VERSION OF THIS!!!! NOT COMPLETE '''                    
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "rest_1", \
                                  "func_motion_correct", \
                                  "rest_calc_tshift_resample_volreg.nii.gz"))
                                    
    # run the workflow
    output = run_func_motion_correct(func_scan, 0, "End", True)

    # make the correlation
    ref_out_data = nb.load(ref_out).get_data()  
    output_data = nb.load(output).get_data()
    
    os.system("rm -R func_motion_correct")


    # create a vector of True and False values
    bool_vector = ref_out_data == output_data

    assert bool_vector.all()



def test_run_functional_brain_mask_3dautomask():

    import os
    import nibabel as nb

    import pkg_resources as p

    from qap.functional_preproc import run_functional_brain_mask


    if "functional_brain_mask" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)

    
    func_motion = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                      "rest_1", \
                                      "func_motion_correct", \
                                      "rest_calc_tshift_resample_" \
                                      "volreg.nii.gz"))
                                    
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "rest_1", \
                                  "functional_brain_mask", \
                                  "rest_calc_tshift_resample_volreg_" \
                                  "mask.nii.gz"))
                                    
    # run the workflow
    output = run_functional_brain_mask(func_motion, False)

    # make the correlation
    ref_out_data = nb.load(ref_out).get_data()  
    output_data = nb.load(output).get_data()
    
    os.system("rm -R functional_brain_mask")


    # create a vector of True and False values
    bool_vector = ref_out_data == output_data

    assert bool_vector.all()



def test_run_functional_brain_mask_BET():

    import os
    import nibabel as nb

    import pkg_resources as p

    from qap.functional_preproc import run_functional_brain_mask


    if "functional_brain_mask" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)

    
    func_motion = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                      "rest_1", \
                                      "func_motion_correct", \
                                      "rest_calc_tshift_resample_" \
                                      "volreg.nii.gz"))
                                    
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "rest_1", \
                                  "functional_brain_mask", \
                                  "rest_calc_tshift_resample_volreg_" \
                                  "mask_BET.nii.gz"))
                                    
    # run the workflow
    output = run_functional_brain_mask(func_motion, True)

    # make the correlation
    ref_out_data = nb.load(ref_out).get_data()  
    output_data = nb.load(output).get_data()
    
    os.system("rm -R functional_brain_mask")


    # create a vector of True and False values
    bool_vector = ref_out_data == output_data

    assert bool_vector.all()



def test_run_mean_functional():

    import os
    import nibabel as nb

    import pkg_resources as p

    from qap.functional_preproc import run_mean_functional


    if "mean_functional" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)

    
    func_motion = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                      "rest_1", \
                                      "func_motion_correct", \
                                      "rest_calc_tshift_resample_" \
                                      "volreg.nii.gz"))
                                    
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "rest_1", \
                                  "mean_functional", \
                                  "rest_calc_tshift_resample_volreg_" \
                                  "tstat.nii.gz"))
                                    
    # run the workflow
    output = run_mean_functional(func_motion)

    # make the correlation
    ref_out_data = nb.load(ref_out).get_data()  
    output_data = nb.load(output).get_data()
    
    os.system("rm -R mean_functional")


    # create a vector of True and False values
    bool_vector = ref_out_data == output_data

    assert bool_vector.all()


    