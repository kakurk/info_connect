% split the whole MTL ROI

if isempty(which('spm'))
    error('Add SPM to the MATLAB searchpath')
end

wholeMTLroiFN = '/gsfs0/data/ritcheym/data/fmri/mice/analysis/retrieval/model_est/lssGenBeta-output/unsmoothed/rMTL_group50_mask.nii';
roipath = '/gsfs0/scratch/kurkela/results/mice-itemret-informational-connectivity/rois';
V = spm_vol(wholeMTLroiFN);
Y = spm_read_vols(V);

%% Left
LEFT   = zeros(size(Y,1), size(Y,2), size(Y,3));
LEFT(1:round(size(LEFT,1)/2),:,:) = Y(1:round(size(LEFT,1)/2),:,:);
LEFT_V = V;
LEFT_V.fname = fullfile(roipath, 'rMTL_group50_mask_L.nii');
spm_write_vol(LEFT_V, LEFT);

%% Right
RIGHT   = zeros(size(Y,1), size(Y,2), size(Y,3));
RIGHT(round(size(LEFT,1)/2):end,:,:) = Y(round(size(LEFT,1)/2):end,:,:);
RIGHT_V = V;
RIGHT_V.fname = fullfile(roipath, 'rMTL_group50_mask_R.nii');
spm_write_vol(RIGHT_V, RIGHT);
