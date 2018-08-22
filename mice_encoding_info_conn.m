function mice_encoding_info_conn(s)

%% Add Paths
% SPM12, CoSMoMVPA, and the Informational Connectivity Toolbox
addpath(genpath('/gsfs0/data/kurkela/Documents/toolboxes-fmri/spm12'));
addpath(genpath('/gsfs0/data/kurkela/Documents/toolboxes-fmri/CoSMoMVPA'));
addpath('/gsfs0/data/kurkela/Documents/latest_IC_toolbox');

%% Relevant Directories
% data_path = where the single trial beta images are
% inital_spm_model_path = where the original SPM models are
% bids_path = where the raw BIDS formatted data are
data_path  = '/gsfs0/data/ritcheym/data/fmri/mice/analysis/encoding/SingleTrialModel_regularmodel';
inital_spm_model_path = '/gsfs0/data/ritcheym/data/fmri/mice/analysis/encoding/Emotion_and_Contexts';

%% Masks
% full brain and MTL masks
maskType = 'HIPP-PHC';
switch maskType
    case 'HIPP-PHC'
        masks{1} = fullfile(data_path, 'rHIPP_BODY_L_mask.nii');
        masks{2} = fullfile(data_path, 'rPHC_ANT_L_mask.nii');
end

%% Subjects
% Subject IDs
%subjects = cellstr(spm_select('List', data_path, 'dir', 'sub-s0[0-3][0-9]'));
subjects  = {'sub-s003' , 'sub-s002', 'sub-s023'};

%% Parameters
% hemodynamic_lag (in TRs).
numberOfROIs        = length(masks);

%% Templates
% In order to perform informational connectivity, we need to define several
% "template patterns" to correlate with the fMRI timecourse. There should
% be one "template pattern" per condition, that is defined based on
% independent training data.

% all of the single trial beta images for this subject
single_trial_betasFNs = RecurseAndFilterFileSearch(data_path, 'Sess0[1-6].*\.nii', subjects{s});

% Grab all four dimensional beta images. Same as single_trial_betasFNs,
% just collapsed into a single file. Four dimensional file format required
% for cosmo_fmri_dataset.m
fourDsingletrialbetasFN = RecurseAndFilterFileSearch(data_path, '.*all-betas\.nii$', subjects{s});
assert(length(fourDsingletrialbetasFN) == 1, 'There is more than 1 four D file');

% initalize
ds_template = cell(1, numberOfROIs);

for m = 1:numberOfROIs
    
    % Load the fMRI data into MATLAB using cosmo_fmri_dataset. Assign a new
    % dataset attribute "maskName" with the name of the mask used to create
    % it.
    ds_template{m} = cosmo_fmri_dataset(char(fourDsingletrialbetasFN), 'mask', masks{m});
    ds_template{m}.a.maskName = masks{m};
    
    % using the custom regularExpression, extract the substring from the
    % single_trial_betaFNs identifying the beta's emotional valence. Take
    % the match and convert it from a nested cell array to a cell string.
    % Assign it as a new sample attribute "EmotionalValence".
    regularExpression = '(?<=Emotion-)[N][a-z]{6,7}(?=_)';
    matches = regexp(single_trial_betasFNs, regularExpression, 'match');
    matches = cellfun(@char, matches, 'UniformOutput', false);
    ds_template{m}.sa.EmotionalValence = matches;
    
    % using the custom regularExpression, extract the substring from the
    % single_trial_betaFNs, identifying the beta's context number. Take the
    % match and convert it from a nested cell array to a cell string to a 
    % vector of doubles. Assign this column vector as a new sample
    % attribute, "ContextNumber".
    regularExpression = '(?<=Context-)[0-4](?=_)';
    matches           = regexp(single_trial_betasFNs, regularExpression, 'match');
    matches           = cellfun(@(x) str2double(cell2mat(x)), matches);
    ds_template{m}.sa.ContextNumber = matches;
    
    % using custom regularExpression, extract the session number from the
    % single_trial_betaFNs.
    regularExpression = '(?<=Sess0)[0-6](?=_)';
    matches           = regexp(single_trial_betasFNs, regularExpression, 'match');
    matches           = cellfun(@(x) str2double(cell2mat(x)), matches);
    ds_template{m}.sa.chunks = matches;
    
end

%% Timecourse
% Use SPM12 and CoSMoMVPA tools to get the fmri data into MATLAB

%%% Timecourse
% Strategy: find the original GLM and figure out the raw timeseries from
% the information contained within the SPM.mat file. Why do it this way? It
% ensures that we are using the exact same timeseries that were used to
% creaete the single trial beta images above.

% Find the original SPM.mat file and load it into MATLAB
SPM_FN = RecurseAndFilterFileSearch(inital_spm_model_path, 'SPM.mat', subjects{s});
assert(length(SPM_FN) == 1, 'More than 1 SPM.mat identified for subject %s', subjects{s})
SPM    = [];
load(char(SPM_FN))

% Figure out the raw timeseries
timeseriesFNs_w_frames = cellstr(SPM.xY.P);
timeseriesFNs   = regexp(timeseriesFNs_w_frames, '.*(?=,[0-9]{1,3})', 'match');
timeseriesFNs   = cellfun(@char, timeseriesFNs, 'UniformOutput', false);
timeseriesFNs   = unique(timeseriesFNs);
timecourse      = cell(1,length(timeseriesFNs));

% for each mask (i.e., ROI)
ds_timecourse = {};
for m = 1:numberOfROIs
    for c = 1:length(timeseriesFNs) % for each run
        
        timecourse{c}            = cosmo_fmri_dataset(timeseriesFNs{c}, ...
                                                      'mask', masks{m}, ...
                                                      'chunks', c);
        timecourse{c}.a.maskName = masks{m};
    end
    ds_timecourse = vertcat(ds_timecourse, cosmo_stack(timecourse));
end

%% Information Connectivity!

keyboard

info_timecourse = cell(1, numberOfROIs);
for m = 1:numberOfROIs
    info_timecourse{m} = cosmo_information_connectivity(ds_template{m}, ds_timecourse{m});
end

keyboard

end