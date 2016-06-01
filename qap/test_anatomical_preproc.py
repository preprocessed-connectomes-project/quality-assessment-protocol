
test_sub_dir = "test_data"


def test_run_anatomical_reorient():

    import os
    import numpy as np
    import nibabel as nb
    import shutil

    import pkg_resources as p

    from qap.anatomical_preproc import run_anatomical_reorient
   
    anatomical_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                          "anatomical_scan.nii.gz"))

    ref_output = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_reorient.nii.gz"))

    out_dir = os.path.join(os.getcwd(), "unit_tests_anat_preproc")
    output = run_anatomical_reorient(anatomical_scan, out_dir=out_dir)

    ref_out_data = nb.load(ref_output).get_data()
    out_data = nb.load(output).get_data()

    try:
        shutil.rmtree(out_dir)
    except:
        pass

    np.testing.assert_array_equal(ref_out_data, out_data)



def test_run_anatomical_skullstrip():

    import os
    import numpy as np
    import nibabel as nb
    import shutil

    import pkg_resources as p

    from qap.anatomical_preproc import run_anatomical_skullstrip
   
    anat_reorient = p.resource_filename("qap", \
                                        os.path.join(test_sub_dir, \
                                        "anat_reorient.nii.gz"))

    ref_output = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_brain.nii.gz"))

    out_dir = os.path.join(os.getcwd(), "unit_tests_anat_preproc")
    output = run_anatomical_skullstrip(anat_reorient, out_dir=out_dir)

    ref_out_data = nb.load(ref_output).get_data()
    out_data = nb.load(output).get_data()

    try:
        shutil.rmtree(out_dir)
    except:
        pass

    np.testing.assert_array_equal(ref_out_data, out_data)



def test_run_afni_anatomical_linear_registration_brain_only():

    import os
    import numpy as np
    import nibabel as nb
    import shutil

    import pkg_resources as p

    from qap.anatomical_preproc import run_afni_anatomical_linear_registration
   
    anat_brain = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_brain.nii.gz"))

    template_brain = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                         "MNI152_T1_3mm_brain.nii.gz"))

    ref_output = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "allineate_warped_brain.nii.gz"))

    out_dir = os.path.join(os.getcwd(), "unit_tests_anat_preproc")
    output = run_afni_anatomical_linear_registration(anat_brain,
                                                     template_brain,
                                                     out_dir=out_dir)

    ref_out_data = nb.load(ref_output).get_data()
    out_data = nb.load(output).get_data()

    try:
        shutil.rmtree(out_dir)
    except:
        pass

    np.testing.assert_array_equal(ref_out_data, out_data)



def test_run_afni_anatomical_linear_registration_skull_on():

    import os
    import numpy as np
    import nibabel as nb
    import shutil

    import pkg_resources as p

    from qap.anatomical_preproc import run_afni_anatomical_linear_registration
   
    anat_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "anat_reorient.nii.gz"))

    template_head = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "MNI152_T1_3mm.nii.gz"))

    ref_output = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "allineate_warped_head.nii.gz"))

    out_dir = os.path.join(os.getcwd(), "unit_tests_anat_preproc")
    output = run_afni_anatomical_linear_registration(anat_reorient,
                                                     template_head,
                                                     skull_on=True,
                                                     out_dir=out_dir)

    ref_out_data = nb.load(ref_output).get_data()
    out_data = nb.load(output).get_data()

    try:
        shutil.rmtree(out_dir)
    except:
        pass

    np.testing.assert_array_equal(ref_out_data, out_data)



def test_run_afni_segmentation():

    import os
    import numpy as np
    import nibabel as nb
    import shutil

    import pkg_resources as p

    from qap.anatomical_preproc import run_afni_segmentation
   
    anat_brain = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_brain.nii.gz"))

    ref_csf_output = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                         "anatomical_csf_mask.nii.gz"))

    ref_gm_output = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "anatomical_gm_mask.nii.gz"))

    ref_wm_output = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "anatomical_wm_mask.nii.gz"))

    out_dir = os.path.join(os.getcwd(), "unit_tests_anat_preproc")
    output = run_afni_segmentation(anat_brain, out_dir=out_dir)

    out_data = {}
    out_ref_data = {}
    for seg_mask in output:
        print seg_mask
        if "csf_mask" in seg_mask:
            out_ref_data["csf"] = nb.load(ref_csf_output).get_data()
            out_data["csf"] = nb.load(seg_mask).get_data()
        elif "gm_mask" in seg_mask:
            out_ref_data["gm"] = nb.load(ref_gm_output).get_data()
            out_data["gm"] = nb.load(seg_mask).get_data()
        elif "wm_mask" in seg_mask:
            out_ref_data["wm"] = nb.load(ref_wm_output).get_data()
            out_data["wm"] = nb.load(seg_mask).get_data()

    try:
        shutil.rmtree(out_dir)
    except:
        pass

    for seg_type in out_data.keys():
        try:
            np.testing.assert_array_equal(out_ref_data[seg_type], \
                                              out_data[seg_type])
        except Exception as e:
            print e, "\n\n", seg_type



def run_all_tests_anatomical_preproc():

    test_run_anatomical_reorient()
    test_run_anatomical_skullstrip()
    test_run_afni_anatomical_linear_registration_brain_only()
    test_run_afni_anatomical_linear_registration_skull_on()
    test_run_afni_segmentation()


