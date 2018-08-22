function IC = mice_encoding_info_conn(s)

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
bids_path  = '/gsfs0/data/ritcheym/data/fmri/mice/data/sourcedata';

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
hemodynamic_lag     = 4;
TR                  = 1.5;
num_of_TRs_per_Sess = 165;

%% Load Data
% Use SPM12 and CoSMoMVPA tools to get the fmri data into MATLAB

%%% Timecourse
% Strategy: find the original GLM and figure out the raw timeseries from
% the information contained within the SPM.mat file.

% Find the original SPM.mat file
dir2search = fullfile(inital_spm_model_path, subjects{s});
SPM_FN   = spm_select('FPList', dir2search, 'SPM.mat');
assert(~isempty(SPM_FN), 'Could not find this subject''s original SPM.mat file', subjects{s})
SPM = [];
load(SPM_FN)

% Figure out the raw timeseries
somefile   = unique(cellfun(@(x) x{:}, regexp(cellstr(SPM.xY.P), '.*(?=,[0-9]{1,3})', 'match'), 'UniformOutput', false));
timecourse = cell(1,length(somefile));

% All of the BIDS information from MICE
BIDS = spm_BIDS(bids_path);

% The BIDS information for this subject
thisSubjectBIDS = BIDS.subjects(strcmp({BIDS.subjects.name}, subjects{s}));

% Just the encoding task
thisSubjectBIDS = thisSubjectBIDS.func(strcmp({thisSubjectBIDS.func.task}, 'encoding'));

% for each mask (i.e., ROI)
concatenated_timecourse = {};
for m = 1:length(masks)
    for c = 1:length(somefile) % for each run
        
        filt{1} = strcmp(thisSubjectBIDS(c).events.ContextNum, '1');
        filt{2} = strcmp(thisSubjectBIDS(c).events.ContextNum, '2');
        filt{3} = strcmp(thisSubjectBIDS(c).events.ContextNum, '3');
        filt{4} = strcmp(thisSubjectBIDS(c).events.ContextNum, '4');
        
        timepointIDs = zeros(num_of_TRs_per_Sess, 1);
        
        for f = 1:length(filt)
            contextNonsetsInTRs = round(thisSubjectBIDS(c).events.onset(filt{f})/TR + hemodynamic_lag);
            timepointIDs(contextNonsetsInTRs)  = f;
        end
        
        timecourse{c}            = cosmo_fmri_dataset(somefile{c}, 'mask', masks{m}, 'chunks', c);
        timecourse{c}.sa.targets = timepointIDs;
        timecourse{c}.sa.roi     = ones(size(timecourse{c}.samples, 1), 1) * m;
        
        %%% Blocks
        blockFilt = strcmp(thisSubjectBIDS(c).events.Condition, 'neg') | ...
                    strcmp(thisSubjectBIDS(c).events.Condition, 'neu');
        % onsets
        blockonsets = thisSubjectBIDS(c).events.onset(blockFilt);
        
        % tricky part: durations. figure out the offset of the fourth trial
        % in each block.
        fourthTrialInEachBlockIDX = find(blockFilt) + 4;
        fourthTrialInEachBlockOnsets = thisSubjectBIDS(c).events.onset(fourthTrialInEachBlockIDX);
        fourthTrialInEachBlockDur    = 3;
        fourthTrialInEachBlockOffset = fourthTrialInEachBlockOnsets + fourthTrialInEachBlockDur;
        blockdurs   = fourthTrialInEachBlockOffset - blockonsets;
        
        timecourse{c}.sa.blockHeatVec = zeros(165, 1);
        
        blockonsetsInTRs = round(blockonsets/TR);
        blockdursInTRs   = round(blockdurs/TR);
        
        blockDurationsAsIDXs = [];
        blockDurationsVIZ = {};
        
        %% A visual
        
        for i = 1:length(blockonsetsInTRs)
            curBlockIDXs = blockonsetsInTRs(i):blockonsetsInTRs(i)+blockdursInTRs(i);
            blockDurationsAsIDXs = horzcat(blockDurationsAsIDXs, curBlockIDXs);
            blockDurationsVIZ = vertcat(blockDurationsVIZ, {curBlockIDXs});
        end
        
        VISUAL = zeros(8, 165);
        
        for i = 1:length(blockDurationsVIZ)
           VISUAL(i, blockDurationsVIZ{i}) = 1; 
        end
        
        %figure; imagesc(VISUAL); colorbar;
        
        %%
        
        timecourse{c}.sa.blockHeatVec(blockDurationsAsIDXs) = 1;

    end
    concatenated_timecourse = vertcat(concatenated_timecourse, cosmo_stack(timecourse));
end

%% Information Connectivity Toolbox
% The information connectivity toolbox requries that the data is structured
% in a very specific way. See run_ROI_IC.m help text for more info.
IC.data       = vertcat(concatenated_timecourse{1}.samples', concatenated_timecourse{2}.samples'); % rows = voxels, columns = trials
A = ones(size(concatenated_timecourse{1}.samples, 2), 1);
B = ones(size(concatenated_timecourse{2}.samples, 2), 1) * 2;
IC.ROIs       = vertcat(A, B);
mat = vertcat(concatenated_timecourse{1}.sa.targets' == 1, ...
              concatenated_timecourse{1}.sa.targets' == 2, ...
              concatenated_timecourse{1}.sa.targets' == 3, ...
              concatenated_timecourse{1}.sa.targets' == 4);
IC.conditions = mat;
IC.folds      = concatenated_timecourse{1}.sa.chunks';
IC.selector   = concatenated_timecourse{1}.sa.blockHeatVec';
IC.ROI_names  = cellfun(@(x) char(regexp(x, 'r[A-Z].*\.nii$', 'match')), masks, 'UniformOutput', false);

IC = run_ROI_IC(IC);

end