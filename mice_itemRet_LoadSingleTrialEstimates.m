function mice_itemRet_LoadSingleTrialEstimates(s)
% load single trial estimates. Write cosmo dataset to disk.

%% Add Paths
% SPM12, CoSMoMVPA, and the Informational Connectivity Toolbox
addpath(genpath('/gsfs0/data/kurkela/Documents/toolboxes-fmri/spm12'));
addpath(genpath('/gsfs0/data/kurkela/Documents/toolboxes-fmri/CoSMoMVPA'));
addpath(genpath('/gsfs0/data/kurkela/Documents/info_conn/thirdparty'));

%% Relevant Directories
% data_path = where the single trial beta images are
% inital_spm_model_path = where the original SPM models are
% bids_path = where the raw BIDS formatted data are
data_path       = '/gsfs0/data/ritcheym/data/fmri/mice/analysis/retrieval/model_est/lssGenBeta-output/unsmoothed';
rosies_roi_path = '/gsfs0/data/ritcheym/data/fmri/mice/analysis/retrieval/model_est/lssGenBeta-output/unsmoothed/ROIs';
roses_roi_path  = '/gsfs0/data/ritcheym/data/fmri/orbit/analysis/orbit/myROIs/MNI-space';
outpath         = '/gsfs0/scratch/kurkela/results/mice-itemret-informational-connectivity';

%% Masks
% full brain and MTL masks

%%% PM-System

    % left
    masks{1}  = fullfile(rosies_roi_path, 'rHIPP_BODY_L_mask.nii');
    masks{2}  = fullfile(rosies_roi_path, 'rPHC_ANT_L_mask.nii');
    masks{3}  = fullfile(roses_roi_path, 'PCC_L_ROI.nii');
    masks{4}  = fullfile(roses_roi_path, 'RSC_L_ROI.nii');
    masks{5}  = fullfile(roses_roi_path, 'PREC_L_ROI.nii');
    masks{6}  = fullfile(roses_roi_path, 'ANG_L_ROI.nii');

    % right
    masks{7}  = fullfile(rosies_roi_path, 'rHIPP_BODY_R_mask.nii');
    masks{8}  = fullfile(rosies_roi_path, 'rPHC_ANT_R_mask.nii');
    masks{9}  = fullfile(roses_roi_path, 'PCC_R_ROI.nii');
    masks{10} = fullfile(roses_roi_path, 'RSC_R_ROI.nii');
    masks{11} = fullfile(roses_roi_path, 'PREC_R_ROI.nii');
    masks{12} = fullfile(roses_roi_path, 'ANG_R_ROI.nii');

%%% AT-System

    % left
    masks{13} = fullfile(rosies_roi_path, 'rAMY_L_mask.nii');
    masks{14} = fullfile(rosies_roi_path, 'rHIPP_HEAD_L_mask.nii');
    masks{15} = fullfile(rosies_roi_path, 'rPRC_L_mask.nii');
    masks{16} = fullfile(roses_roi_path, 'ITC_L_ROI.nii');
    masks{17} = fullfile(roses_roi_path, 'OFC_L_ROI.nii');

    % right
    masks{18} = fullfile(rosies_roi_path, 'rAMY_R_mask.nii');
    masks{19} = fullfile(rosies_roi_path, 'rHIPP_HEAD_R_mask.nii');
    masks{20} = fullfile(rosies_roi_path, 'rPRC_R_mask.nii');
    masks{21} = fullfile(roses_roi_path, 'ITC_R_ROI.nii');
    masks{22} = fullfile(roses_roi_path, 'OFC_R_ROI.nii');

%%% MTL
    
    masks{23}  = fullfile(outpath, 'rois', 'rMTL_group50_mask_L.nii');
    masks{24}  = fullfile(outpath, 'rois', 'rMTL_group50_mask_R.nii');


outpath = fullfile(outpath, 'SingleTrialEstimates');
if ~exist(outpath, 'dir')
    mkdir(outpath)
end

%% Subjects
% Subject IDs
subjects = cellstr(spm_select('List', data_path, 'dir', 's0[0-3][0-9]'));
numberOfROIs = length(masks);

%% Load Single Trial Estimates

% all of the single trial beta images for this subject
single_trial_betasFNs = RecurseAndFilterFileSearch(data_path, 'Sess00[1-6].*\.nii', [filesep subjects{s}]);
single_trial_betasFNs(~contains(single_trial_betasFNs, 'itemret')) = [];
single_trial_betasFNs(contains(single_trial_betasFNs, '-OLD')) = [];

% We need to 4-d file them.
fourDsingletrialbetasFN = threeDTofourD(single_trial_betasFNs, fileparts(outpath));

% initalize
ds = cell(1, numberOfROIs);

for m = 1:numberOfROIs
    
    % see cosmo_fmri_dataset
    ds{m} = cosmo_fmri_dataset(char(fourDsingletrialbetasFN), 'mask', masks{m});
    
    % assign a new feature attribute (.fa) to label each voxel in this
    % dataset as belonging to this ROI
    
    % derive the ROI label from the mask filename
    [~, ROIlabel, ~] = fileparts(masks{m});
    
    % remove the leading 'r' and the trailing '_mask' from the ROIlabel
    ROIlabel = regexprep(ROIlabel, '^r', '');
    ROIlabel = regexprep(ROIlabel, '_mask', '');
    
    % feature attribute must be the same size as the features
    ds{m}.fa.ROIlabel = repmat({ROIlabel}, 1, size(ds{m}.samples, 2));
    
end

% see cosmo_stack
ds = cosmo_stack(ds, 2);

% Assign Dataset Attributes
ds.a.ExperimentName = 'MICE_itemret';
ds.a.SubjectID      = subjects{s};

% using the custom regularExpression, extract the substring from the
% single_trial_betaFNs identifying the beta's Probe status (i.e., 
% Target/Lure/Novel). Take the match and convert it from a nested cell 
% array to a cell string. Assign it as a new sample attribute "Probe".
regularExpression    = '[tln][auo][rv][ge][l]?(?=-)';
matches              = regexp(single_trial_betasFNs, regularExpression, 'match');
matches              = cellfun(@char, matches, 'UniformOutput', false);
ds.sa.Probe = matches;

% using the custom regularExpression, extract the substring from the
% single_trial_betaFNs identifying the beta's Memory status (i.e., 
% Hit/Miss/FA/CR). Take the match and convert it from a nested cell 
% array to a cell string. Assign it as a new sample attribute 
% "Memory".
regularExpression     = '(?<=-)[hmCF][iRA][ts]?[s]?(?=_)';
matches               = regexp(single_trial_betasFNs, regularExpression, 'match');
matches               = cellfun(@char, matches, 'UniformOutput', false);
ds.sa.Memory = matches;

% using custom regularExpression, extract the session number from the
% single_trial_betaFNs. Take the match and convert it from a nested cell
% array to a double vector. Assign it to a new sample attribute "chunks".
regularExpression     = '(?<=Sess00)[0-6](?=_)';
matches               = regexp(single_trial_betasFNs, regularExpression, 'match');
matches               = cellfun(@(x) str2double(cell2mat(x)), matches);
ds.sa.chunks = matches;

% Now for the context information...which is not contained within the
% beta file names :(. Strategy: find the singleTrial information *.mat
% files created by Maureen's generate_single_trial.m script. Load it up and
% extract out the context information and assign it as a sample attribute
% "ContextNum".
singleTrialInfoMat = RecurseAndFilterFileSearch(data_path, 'MICE_singletrial.*\.mat', subjects{s});
load(char(singleTrialInfoMat), 'subTable');
function out = custom_extraction(x)
    if ~isempty(x)
        out = cell2mat(x);
    else
        out = NaN;
    end
end
ds.sa.ContextNum = cellfun(@custom_extraction, subTable.context(contains(subTable.trial, 'itemret')));

% Like above with the context information, extract the emotional valence of
% the current trial from the "codnition" field of the singleTrial
% information table created by Maureen's generate_single_trial.m script.
% Extract the information and assign it as a sample attribute
% "EmotionalValence"
function out = custom_extraction_emo(x)
    if ismember(x, [1 2])
        out = {'neg'};
    elseif ismember(x, [3 4])
        out = {'neut'};
    end
end
ds.sa.EmotionalValence = cellfun(@custom_extraction_emo, subTable.condition(contains(subTable.trial, 'itemret')));

outfilename = fullfile(outpath, sprintf('sub-%s_SingleTrialEstimates.mat', subjects{s}));
save(outfilename, 'ds')

%% Subfunctions 

function fullpathtofourDfile = threeDTofourD(files, outpath)

    fourD_filename = sprintf('sub-%s_itemret_betas.nii', subjects{s});
    
    fullpathtofourDfile = fullfile(outpath, fourD_filename);
    
    if ~exist(char(fullpathtofourDfile), 'file')
    
        matlabbatch{1}.spm.util.cat.vols  = files;
        matlabbatch{1}.spm.util.cat.name  = fullfile(outpath, fourD_filename);
        matlabbatch{1}.spm.util.cat.dtype = 0;

        spm_jobman('run', matlabbatch)
    
    end
    
end
end