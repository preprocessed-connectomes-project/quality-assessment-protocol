import nibabel as nib
import numpy as np
from scipy import stats


def remove_zero_variance_voxels(func_timeseries, mask):

    for i in range(0, len(func_timeseries)):

        for j in range(0, len(func_timeseries[0])):

            for k in range(0, len(func_timeseries[0][0])):

                var = func_timeseries[i][j][k].var()

                if int(var) == 0:

                    mask[i][j][k] = mask[i][j][k] * 0

    return mask



def load(func_file, mask_file, check4d=True):

    func_img    = nib.load(func_file)
    mask_img    = nib.load(mask_file)

    mask        = mask_img.get_data()
    func        = func_img.get_data().astype(np.float)

    if check4d and len(func.shape) != 4:
        raise Exception("Input functional %s should be 4-dimensional" % func_file)

    mask_var_filtered = remove_zero_variance_voxels(func, mask)

    func        = func[mask_var_filtered.nonzero()].T # will have ntpts x nvoxs
    
    
    return func



def robust_stdev(func, interp="fraction"):
    """
    Compute robust estimation of standard deviation
    """
    lower_qs    = np.percentile(func, 25, axis=0)
    upper_qs    = np.percentile(func, 75, axis=0)
    # note: won't work on roxy with scipy == 0.9
    #lower_qs    = stats.scoreatpercentile(func, 25, interpolation_method=interp, axis=0)
    #upper_qs    = stats.scoreatpercentile(func, 75, interpolation_method=interp, axis=0)
    stdev       = (upper_qs - lower_qs)/1.349
    return stdev



def ar_nitime(x, order=1, center=False):
    """
    Borrowed from nipy.algorithms.AR_est_YW.
    aka from nitime import algorithms as alg.
    
    We could speed this up by having the autocorr only compute lag1.
    """
    from nitime.lazy import scipy_linalg as linalg
    import nitime.utils as utils
    if center:
        x = x.copy()
        x = x - x.mean()
    r_m = utils.autocorr(x)[:order + 1]
    Tm  = linalg.toeplitz(r_m[:order])
    y   = r_m[1:]
    ak  = linalg.solve(Tm, y)
    return ak[0]



def ar_statsmodels(x, order=(1,0), trend='nc'):
    import statsmodels.api as sm
    arma_mod = sm.tsa.ARMA(x)
    arma_res = arma_mod.fit(order=order, trend=trend, disp=False)
    return arma_res.arparams[0]



def ar1(func, method=ar_nitime):
    func_centered = func - func.mean(0)
    #import code
    #code.interact(local=locals())
    ar_vals = np.apply_along_axis(method, 0, func_centered)
    return ar_vals



def calc_dvars(func, output_all=False, interp="fraction"):
    # Robust standard deviation
    func_sd     = robust_stdev(func, interp)
    
    # AR1
    func_ar1    = ar1(func)

    # Predicted standard deviation of temporal derivative
    func_sd_pd  = np.sqrt(2 * (1 - func_ar1)) * func_sd
    diff_sd_mean= func_sd_pd.mean()

    # Compute temporal difference time series 
    func_deriv  = np.diff(func, axis=0)

    # DVARS
    ## (no standardization)
    dvars_plain = func_deriv.std(1, ddof=1) # TODO: Why are we not ^2 this & getting the sqrt?
    ## standardization
    dvars_stdz  = dvars_plain/diff_sd_mean
    ## voxelwise standardization
    diff_vx_stdz= func_deriv/func_sd_pd
    dvars_vx_stdz = diff_vx_stdz.std(1, ddof=1)
    
    if output_all:
        out = np.vstack((dvars_stdz, dvars_plain, dvars_vx_stdz))
    else:
        out = dvars_stdz.reshape(len(dvars_stdz), 1)
    
    return out



def calc_mean_dvars(dvars):
    mean_dvars = dvars.mean(0)
    return mean_dvars



def mean_dvars_wrapper(func_file, mask_file, dvars_out_file=None):
    func    = load(func_file, mask_file)
    dvars   = calc_dvars(func)
    if dvars_out_file:
        np.savetxt(dvars_out_file, dvars, fmt='%.12f')
    mean_d  = calc_mean_dvars(dvars)
    return mean_d[0]



def test():
    func    = load("sample_func.nii.gz", "sample_func_mask.nii.gz")
    dvars   = calc_dvars(func)
    mean_d  = calc_mean_dvars(dvars)
    ref_dvars = np.loadtxt("sample_dvars.txt")
    ref_dvars = ref_dvars[:,[1,0,2]]
    


def specific_tests():
    import numpy.testing
    
    ffile   = "sample_func.nii.gz"
    mfile   = "sample_func_mask.nii.gz"
    func    = load(ffile, mfile)
    
    # Robust Standard Deviation
    ## Differences in the two approaches exist because python will handle
    ## ties via averaging
    func_sd     = robust_stdev(func)
    ref_sd      = load("ref_dvars/DVARS-1033--SD.nii.gz", mfile, False)
    print np.abs(func_sd - ref_sd).mean()
    ## so we can fix this issue by changing the interpolation/ties approach
    ## note however that normally we will use interp="fraction"
    func_sd     = robust_stdev(func, interp="higher")
    ## test
    assert_allclose(func_sd, ref_sd)
    print np.abs(func_sd - ref_sd).mean()
    
    # AR1
    func_ar1    = ar1(func)
    ref_ar1     = load("ref_dvars/DVARS-1033--AR1.nii.gz", mfile, False)
    ## test
    assert_allclose(func_ar1, ref_ar1, 1e-6, 1e-6)
    print np.abs(func_ar1 - ref_ar1).mean()
    
    # Predicted standard deviation of temporal derivative
    func_sd_pd  = np.sqrt(2 * (1 - func_ar1)) * func_sd
    diff_sd_mean= func_sd_pd.mean()
    ref_sd_pd   = load("ref_dvars/DVARS-1033--DiffSDhat.nii.gz", mfile, False)
    ## test
    assert_allclose(func_sd_pd, ref_sd_pd, 1e-5, 1e-5)
    print np.abs(func_sd_pd - ref_sd_pd).mean()
    
    # Compute temporal difference time series
    func_deriv  = np.diff(func, axis=0)
    ref_deriv   = load("ref_dvars/DVARS-1033--Diff.nii.gz", mfile)
    ## test (these should be flipped in sign)
    assert_equal(func_deriv, -1*ref_deriv)
    print np.abs(func_deriv - (-1*ref_deriv)).mean()
    
    # DVARS (no standardization)
    dvars_plain = func_deriv.std(1, ddof=1)
    ref_plain   = np.loadtxt("ref_dvars/DVARS-1033--DiffSD.dat")
    print np.abs(dvars_plain - ref_plain).mean()
    print np.vstack((dvars_plain[:10], ref_plain[:10])).T
    print np.vstack((func_deriv.std(1, ddof=1)[:10], ref_plain[:10])).T
    ## it seems like the differences above stem from an incorrect averaging by my modication of DVARS
    func_img    = nib.load("ref_dvars/DVARS-1033--Diff.nii.gz")
    func_all    = func_img.get_data().astype(np.float)
    func_all    = func_all.reshape(np.prod(func_all.shape[:3]), func_all.shape[-1])
    func_all    = func_all.T
    print func_all[0,func_all[0,:].nonzero()].std(ddof=1) - ref_plain[0]
    
    # DVARS
    ## (no standardization)
    dvars_plain = func_deriv.std(1, ddof=1) # TODO: Why are we not ^2 this & getting the sqrt?
    ## standardization
    dvars_stdz  = dvars_plain/diff_sd_mean
    ## voxelwise standardization
    diff_vx_stdz= func_deriv/func_sd_pd
    dvars_vx_stdz = diff_vx_stdz.std(1, ddof=1)
    ## reference
    ref_dvars = np.loadtxt("sample_dvars.txt")
    ref_dvars = ref_dvars[:,[1,0,2]]
    ## test
    print np.abs(dvars_plain - ref_dvars[:,0]).mean()
    
## Below tests show that ar_nitime is much faster
## so we will use that
# %timeit ar_nitime(func_centered[:,0])
# %timeit ar_statsmodels(func_centered[:,0])
