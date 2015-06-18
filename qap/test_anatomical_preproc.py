

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



def run_all_tests_anatomical_preproc():

    test_run_anatomical_reorient()
    test_run_anatomical_skullstrip()
    
    

