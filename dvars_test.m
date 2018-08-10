addpath("../DVARS");

func_reorient = "qap/test_data/func_reorient.nii.gz";
func_mask = "qap/test_data/fsl_bet_mask.nii.gz";
tol=1e-6;

% read and conform the functional data
Vfunc = niftiread(func_reorient);
Vfunc_sz = size(Vfunc);
Vfunc = double(reshape(Vfunc, prod(Vfunc_sz(1:3)), Vfunc_sz(4)));

% read and conform the mask
Vmask = niftiread(func_mask);
Vmask_sz = size(Vmask);
Vmask = reshape(Vmask, prod(Vmask_sz(1:3)), 1);

% binarize the mask
Vmask(Vmask ~= 0) = 1;
Vmask(Vmask ~= 1) = 0;

% make sure there are no voxels with zero variance
Vmask(abs(var(Vfunc, 0, 2)) < tol) = 0;

% reduce functional data to just the voxels in mask
Vfunc = Vfunc(Vmask == 1, :);

% calculate the DVARS using the code from github, calc scale so mean is 100
[DVARS,DVARS_Stat]=DVARSCalc(Vfunc,'scale',100/mean(Vfunc(:)),'TransPower',1/3,'RDVARS','verbose',1);

%
fid = fopen('/tmp/dvars_matlab_result.json','w');
fprintf(fid,jsonencode(DVARS_Stat));
fclose(fid);
