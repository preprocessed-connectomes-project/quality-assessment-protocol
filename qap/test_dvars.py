
import pytest
import unittest
test_sub_dir = "test_data"


class TestStandardizedDVARS(unittest.TestCase):

    def setUp(self):
        import os
        import pkg_resources as p

        self.func_reorient = \
            p.resource_filename("qap", os.path.join(test_sub_dir,
                                                    "func_reorient.nii.gz"))
        self.func_mask = \
            p.resource_filename("qap",
                                os.path.join(test_sub_dir,
                                             "functional_brain_mask.nii.gz"))

        # this is based on the 2/19/2017 post-fix version of DVARS.sh with
        # standard deviation replaced with root mean square
        self.ref_dvars = [[1.19471], [1.74488], [0.973054], [0.998598],
                          [1.74141], [1.14707], [1.06522], [1.75126],
                          [0.821287], [1.13746], [3.34425], [1.0205],
                          [2.59632], [0.88072], [0.652472], [0.653099],
                          [0.753194], [0.777735], [0.787621], [0.803559],
                          [0.624143], [0.568267], [0.591], [0.599311],
                          [0.673767], [0.707941], [0.716991], [0.791894],
                          [0.748472], [1.07407], [2.67802], [0.816987],
                          [2.45072], [0.763636], [0.760437], [0.884735],
                          [0.794975], [0.966546], [0.670729], [0.685647],
                          [1.02114], [0.769986], [0.90958], [0.870556],
                          [0.995202], [0.801353], [0.976579], [0.825581],
                          [0.708622], [0.69978], [0.58894], [0.588187],
                          [0.772176], [0.723978], [0.717768], [0.825696],
                          [0.741694], [0.61278], [0.685251], [0.639529],
                          [0.695137], [0.67272], [0.814177], [0.777172],
                          [1.07881], [1.07745], [2.71734], [1.02728],
                          [2.4903], [0.832773], [0.719], [0.8396], [0.658567],
                          [0.58555], [0.566178], [0.59842], [0.699388],
                          [0.743473], [0.761807], [0.76259], [3.38393],
                          [0.937676], [2.95067], [1.19944], [1.06869],
                          [1.14729], [0.849156], [0.709569], [0.693313],
                          [0.698364], [0.693688], [0.670661], [0.788288],
                          [0.783307], [0.59745], [0.653137], [0.890101],
                          [1.345], [0.702011], [0.820606], [0.703725],
                          [0.677232], [0.576504], [0.589431], [0.60131],
                          [0.648679], [0.709758], [1.05289], [0.964687],
                          [0.704595], [1.82044], [0.814146], [1.20721],
                          [0.680779], [0.922978], [1.06627], [0.974431],
                          [0.800713], [0.783288], [0.630774], [0.623251],
                          [1.26908], [0.746998]]

        '''
        self.ref_dvars = [[1.61668], [1.55594], [1.50794], [1.49141],
                          [2.2812], [1.69892], [1.70481], [1.76283],
                          [1.02876], [1.58089], [1.10771], [1.28999],
                          [1.32713], [1.01431], [0.842966], [0.943889],
                          [1.27805], [1.13388], [1.23735], [1.23197],
                          [1.00057], [0.941476], [0.951362], [0.970872],
                          [1.16752], [1.06146], [1.06044], [0.977928],
                          [1.0234], [1.12204], [0.938907], [0.845179],
                          [0.941355], [1.04987], [1.09318], [0.978873],
                          [1.20858], [1.29815], [1.05929], [0.990305],
                          [1.25383], [1.063], [1.11785], [1.0264], [0.961176],
                          [1.14121], [1.77877], [1.39043], [1.20306],
                          [1.24328], [0.975394], [0.956231], [1.10099],
                          [1.13488], [1.17661], [1.49438], [1.30942],
                          [0.906392], [0.964955], [0.978622], [1.02638],
                          [0.976282], [1.2036], [0.991174], [1.3067],
                          [1.22121], [1.30675], [1.19907], [1.24118],
                          [1.0608], [1.11127], [1.08457], [0.97589],
                          [0.953487], [0.919849], [1.01331], [1.06032],
                          [1.04153], [0.945779], [0.919424], [0.881916],
                          [0.850963], [0.97861], [0.966441], [1.17289],
                          [1.24443], [1.05687], [1.06819], [1.10116],
                          [1.02466], [0.953319], [0.922557], [0.950052],
                          [0.904365], [0.976654], [1.05072], [1.09395],
                          [1.0583], [1.0448], [1.08735], [1.01156], [1.11197],
                          [0.891858], [0.977091], [0.962705], [1.00625],
                          [1.01554], [1.03546], [1.04511], [1.09378],
                          [1.18882], [1.05104], [0.981128], [0.916849],
                          [1.16562], [1.10247], [1.17634], [1.04704],
                          [1.09401], [1.00456], [1.05005], [1.2189],
                          [1.23903]]
        '''

    def test_calc_dvars(self):
        import numpy.testing as nt
        from qap.dvars import calc_dvars
        dvars_out = calc_dvars(self.func_reorient, self.func_mask)
        nt.assert_array_almost_equal(dvars_out, self.ref_dvars)


@pytest.mark.skip
@pytest.mark.quick
def test_remove_zero_variance_voxels():

    import os
    import pickle
    import pkg_resources as p
    
    import nibabel as nb
    import numpy as np
    
    from qap.dvars import remove_zero_variance_voxels

    func_reorient = p.resource_filename("qap", os.path.join(test_sub_dir,
                                        "func_reorient.nii.gz"))
                                      
    func_mask = p.resource_filename("qap", os.path.join(test_sub_dir,
                                    "functional_brain_mask.nii.gz"))
                                    
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir,
                                  "no_zero_variance_voxels_mask.p"))
                                   
    func_img = nb.load(func_reorient)
    mask_img = nb.load(func_mask)
    
    func_data = func_img.get_data()
    mask_data = mask_img.get_data()
                                    
    out_mask_data = remove_zero_variance_voxels(func_data, mask_data)
                                    
    with open(ref_out, "r") as f:
        ref_mask_data = pickle.load(f)
        
    np.testing.assert_array_equal(ref_mask_data, out_mask_data)


@pytest.mark.skip
@pytest.mark.quick
def test_load():

    import os
    import pickle
    import pkg_resources as p
    import numpy as np
    
    from qap.dvars import load

    func_reorient = p.resource_filename("qap", os.path.join(test_sub_dir,
                                        "func_reorient.nii.gz"))
                                      
    func_mask = p.resource_filename("qap", os.path.join(test_sub_dir,
                                    "functional_brain_mask.nii.gz"))
                                    
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir,
                                  "loaded_func.p"))
                                    
    func_out_data = load(func_reorient, func_mask)

    # match the reference array
    func_out_data = func_out_data[0:10]
   
    with open(ref_out, "r") as f:
        ref_out_data = pickle.load(f)

    np.testing.assert_array_equal(ref_out_data, func_out_data)
    

@pytest.mark.skip
@pytest.mark.quick
def test_robust_stdev():

    import os
    import pickle
    import pkg_resources as p
    import numpy as np
    
    from qap.dvars import robust_stdev

    func_data_file = p.resource_filename("qap", os.path.join(test_sub_dir,
                                         "loaded_func.p"))
                                                                         
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir,
                                  "robust_stdev_output.p"))
                                    
    with open(func_data_file, "r") as f:
        func_data = pickle.load(f)
    
    with open(ref_out, "r") as f:
        ref_mask_data = pickle.load(f)        
    
    func_out_data = robust_stdev(func_data)
        
    np.testing.assert_array_equal(ref_mask_data, func_out_data)


@pytest.mark.skip
@pytest.mark.quick
def test_ar1():

    import os
    import pickle
    import pkg_resources as p
    import numpy as np
    
    from qap.dvars import ar1

    func_data_file = p.resource_filename("qap", os.path.join(test_sub_dir,
                                                             "loaded_func.p"))
                                                                         
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir,
                                                      "ar1_output.p"))
                                    
    with open(func_data_file, "r") as f:
        func_data = pickle.load(f)
    
    with open(ref_out, "r") as f:
        ref_out_data = pickle.load(f)
        
    func_out_data = ar1(func_data)
        
    np.testing.assert_array_almost_equal(ref_out_data, func_out_data)
