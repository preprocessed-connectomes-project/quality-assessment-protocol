
import pytest
test_sub_dir = "test_data"

@pytest.mark.quick
def test_run_everything_qap_anatomical_spatial_workflow_graph():

    # this tests the workflow builder, not the end results (these have their
    # own unit tests)

    import os
    import shutil

    import pkg_resources as p

    from qap.qap_workflows import run_everything_qap_anatomical_spatial

    anatomical_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                          "anatomical_scan.nii.gz"))

    template_head = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "MNI152_T1_3mm.nii.gz"))

    ref_graph = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anatomical_spatial_graph.dot"))

    out_dir = os.path.join(os.getcwd(), "unit_tests_qap_workflows")
    out_workflow = run_everything_qap_anatomical_spatial(anatomical_scan, \
                                                         template_head, \
                                                         "participant_1",
                                                         out_dir=out_dir, \
                                                         run=False)

    out_workflow_obj = out_workflow[0]
    out_workflow_dir = out_workflow[1]

    # write the dependency graph of the workflow we are testing
    out_graph = os.path.join(out_workflow_dir, "graph.dot")
    out_workflow_obj.write_graph(dotfilename=out_graph, simple_form=False)

    # load the both the reference and the to-test dependency graphs
    with open(ref_graph,"r") as f:
        ref_graph_lines = sorted(f.readlines())

    with open(out_graph,"r") as f:
        out_graph_lines = sorted(f.readlines())

    try:
        shutil.rmtree(out_dir)
    except:
        pass

    assert ref_graph_lines == out_graph_lines


@pytest.mark.quick
def test_run_everything_qap_functional_spatial_workflow_graph():

    # this tests the workflow builder, not the end results (these have their
    # own unit tests)

    import os
    import shutil

    import pkg_resources as p

    from qap.qap_workflows import run_everything_qap_functional_spatial

    functional_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                          "functional_scan.nii.gz"))

    ref_graph = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "functional_spatial_graph.dot"))

    out_dir = os.path.join(os.getcwd(), "unit_tests_qap_workflows")
    out_workflow = run_everything_qap_functional_spatial(functional_scan, \
                                                         "participant_1",
                                                         out_dir=out_dir, \
                                                         run=False)

    out_workflow_obj = out_workflow[0]
    out_workflow_dir = out_workflow[1]

    # write the dependency graph of the workflow we are testing
    out_graph = os.path.join(out_workflow_dir, "graph.dot")
    out_workflow_obj.write_graph(dotfilename=out_graph, simple_form=False)

    # load the both the reference and the to-test dependency graphs
    with open(ref_graph,"r") as f:
        ref_graph_lines = sorted(f.readlines())

    with open(out_graph,"r") as f:
        out_graph_lines = sorted(f.readlines())

    try:
        shutil.rmtree(out_dir)
    except:
        pass

    assert ref_graph_lines == out_graph_lines


@pytest.mark.quick
def test_run_everything_qap_functional_temporal_workflow_graph():

    # this tests the workflow builder, not the end results (these have their
    # own unit tests)

    import os
    import numpy as np
    import nibabel as nb
    import shutil

    import pkg_resources as p

    from qap.qap_workflows import run_everything_qap_functional_temporal

    functional_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                          "functional_scan.nii.gz"))

    ref_graph = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "functional_temporal_graph.dot"))

    out_dir = os.path.join(os.getcwd(), "unit_tests_qap_workflows")
    out_workflow = run_everything_qap_functional_temporal(functional_scan, \
                                                          "participant_1",
                                                          out_dir=out_dir, \
                                                          run=False)

    out_workflow_obj = out_workflow[0]
    out_workflow_dir = out_workflow[1]

    # write the dependency graph of the workflow we are testing
    out_graph = os.path.join(out_workflow_dir, "graph.dot")
    out_workflow_obj.write_graph(dotfilename=out_graph, simple_form=False)


    # load the both the reference and the to-test dependency graphs
    with open(ref_graph,"r") as f:
        ref_graph_lines = sorted(f.readlines())

    with open(out_graph,"r") as f:
        out_graph_lines = sorted(f.readlines())

    try:
        shutil.rmtree(out_dir)
    except:
        pass

    assert ref_graph_lines == out_graph_lines
