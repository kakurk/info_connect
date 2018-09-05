function mice_encoding_info_conn(s)

%% Add Paths
% SPM12, CoSMoMVPA, and the Informational Connectivity Toolbox
addpath(genpath('/gsfs0/data/kurkela/Documents/toolboxes-fmri/spm12'));
addpath(genpath('/gsfs0/data/kurkela/Documents/toolboxes-fmri/CoSMoMVPA'));
addpath(genpath('/gsfs0/data/kurkela/Documents/info_conn/thirdparty'));

%% Relevant Directories
% data_path = where the single trial beta images are
% inital_spm_model_path = where the original SPM models are
% bids_path = where the raw BIDS formatted data are
data_path  = '/gsfs0/data/ritcheym/data/fmri/mice/analysis/encoding/SingleTrialModel_regularmodel';
denoisedRawDataPath = '/gsfs0/data/ritcheym/data/fmri/mice/analysis/encoding/Denoised_Raw_Data';
behav_data_path = '/gsfs0/data/ritcheym/data/fmri/mice/data/sourcedata';
outpath = '/gsfs0/scratch/kurkela/results/mice-encoding-information-connectivity';

%% Masks
% full brain and MTL masks
maskType = 'PM-System';
switch maskType
    case 'PM-System'
        % left
        masks{1} = fullfile(data_path, 'rHIPP_BODY_L_mask.nii');
        masks{2} = fullfile(data_path, 'rPHC_ANT_L_mask.nii');
        % right
        masks{3} = fullfile(data_path, 'rHIPP_BODY_R_mask.nii'); 
        masks{4} = fullfile(data_path, 'rPHC_ANT_R_mask.nii');
    case 'AT-System'
        % left
        masks{1} = fullfile(data_path, 'rAMY_L_mask.nii');
        masks{2} = fullfile(data_path, 'rHIPP_HEAD_L_mask.nii');
        masks{3} = fullfile(data_path, 'rPRC_L_mask.nii');
        % right
        masks{4} = fullfile(data_path, 'rAMY_R_mask.nii');
        masks{5} = fullfile(data_path, 'rHIPP_HEAD_R_mask.nii');
        masks{6} = fullfile(data_path, 'rPRC_R_mask.nii');
end

%% Subjects
% Subject IDs
subjects = cellstr(spm_select('List', data_path, 'dir', 'sub-s0[0-3][0-9]'));
%subjects  = {'sub-s003' , 'sub-s002', 'sub-s023'};

%% Parameters
% hemodynamic_lag (in TRs).
numberOfROIs        = length(masks);
TR                  = 1.5;
TRsPerSession       = 165;
information = 'EmotionalValence'; % 'ContextNumber'

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
    
    % see cosmo_fmri_dataset
    ds_template{m} = cosmo_fmri_dataset(char(fourDsingletrialbetasFN), 'mask', masks{m});
    
    % assign a new feature attribute (.fa) to label each voxel in this
    % dataset as belonging to this ROI
    
    % derive the ROI label from the mask filename
    [~, ROIlabel, ~] = fileparts(masks{m});
    
    % remove the leading 'r' and the trailing '_mask' from the ROIlabel
    ROIlabel = regexprep(ROIlabel, '^r', '');
    ROIlabel = regexprep(ROIlabel, '_mask', '');
    
    % feature attribute must be the same size as the features
    ds_template{m}.fa.ROIlabel = repmat({ROIlabel}, 1, size(ds_template{m}.samples, 2));
    
end

% see cosmo_stack
ds_template = cosmo_stack(ds_template, 2);

% Assign Dataset Attributes
ds_template.a.ExperimentName = 'MICE_encoding';
ds_template.a.SubjectID      = subjects{s};

% using the custom regularExpression, extract the substring from the
% single_trial_betaFNs identifying the beta's emotional valence. Take
% the match and convert it from a nested cell array to a cell string.
% Assign it as a new sample attribute "EmotionalValence".
regularExpression = '(?<=Emotion-)[N][a-z]{6,7}(?=_)';
matches = regexp(single_trial_betasFNs, regularExpression, 'match');
matches = cellfun(@char, matches, 'UniformOutput', false);
ds_template.sa.EmotionalValence = matches;

% using the custom regularExpression, extract the substring from the
% single_trial_betaFNs, identifying the beta's context number. Take the
% match and convert it from a nested cell array to a cell string to a 
% vector of doubles. Assign this column vector as a new sample
% attribute, "ContextNumber".
regularExpression = '(?<=Context-)[0-4](?=_)';
matches           = regexp(single_trial_betasFNs, regularExpression, 'match');
matches           = cellfun(@(x) str2double(cell2mat(x)), matches);
ds_template.sa.ContextNumber = matches;

% using custom regularExpression, extract the session number from the
% single_trial_betaFNs.
regularExpression = '(?<=Sess0)[0-6](?=_)';
matches           = regexp(single_trial_betasFNs, regularExpression, 'match');
matches           = cellfun(@(x) str2double(cell2mat(x)), matches);
ds_template.sa.chunks = matches;

%% Timecourses
% Use SPM12 and CoSMoMVPA tools to get the fmri data into MATLAB

% Residual Images are what is "left over" after accounting for motion and
% session means. SPM mat as information about the original raw data.
ResidualImages = RecurseAndFilterFileSearch(denoisedRawDataPath, '^Res_[0-9]{4}.*\.nii', subjects{s});
SPM_mat_FN     = RecurseAndFilterFileSearch(denoisedRawDataPath, 'SPM.mat', subjects{s});

% Get the full paths to the original raw data from the SPM.mat
SPM = [];
load(char(SPM_mat_FN))
rawDataFNs = cellstr(SPM.xY.P);

% ResidualImages 3-D --> 4-D. They need to be a single 4-D .nii file
ResuidualImages4D = threeDTofourD(ResidualImages);

% grab denoised timecourses for each ROI
for m = 1:numberOfROIs

    % see cosmo_fmri_dataset
    ds_timecourse{m} = cosmo_fmri_dataset(char(ResuidualImages4D), 'mask', masks{m});

    % Add ROIlabel feature attribute
    [~, ROIlabel, ~] = fileparts(masks{m});
    
    % remove the leading 'r' and the trailing '_mask'
    ROIlabel = regexprep(ROIlabel, '^r', '');
    ROIlabel = regexprep(ROIlabel, '_mask', '');
    
    % feature attributes must be the same length as the features
    ds_timecourse{m}.fa.ROIlabel = repmat({ROIlabel}, 1, size(ds_timecourse{m}.samples, 2));
        
end

% Stack the features
ds_timecourse = cosmo_stack(ds_timecourse, 2);

%%% Assign Dataset Attributes
ds_timecourse.a.ExperimentName = 'MICE_encoding';
ds_timecourse.a.SubjectID      = subjects{s};

%%% Assign Sample Attributes

% chunks = runs = scanning sessions
chunks = cellfun(@(x) str2double(char(x)), regexp(rawDataFNs, '(?<=run-)0[0-6]', 'match'));
ds_timecourse.sa.chunks = chunks;

% z-score with each scanning session
ds_timecourse = cosmo_fx(ds_timecourse, @zscore, {'chunks'});

%---Read in Behavioral Data -----%

% identifty the events tsv files
behav_data = RecurseAndFilterFileSearch(behav_data_path, '.*encoding.*events\.tsv', subjects{s});

% create a logical filter, identifying all of the runs that were included
% in the model
chunk_filt = false(length(behav_data), 1);
for c = unique(chunks)'
    filt = contains(behav_data, sprintf('run-0%d', c));
    chunk_filt = filt | chunk_filt;
end

% censor runs that were not in the model
behav_data(~chunk_filt) = [];

% read into matlab
behav_data = cellfun(@custom_read_tsv, behav_data, 'UniformOutput', false);
behav_data = vertcat(behav_data{:});

% EmotionalValence = the emotional valence that is occuring at each
% timepoint. n/a when fixation is on the screen.
ds_timecourse.sa.EmotionalValence = repmat({'n/a'}, size(ds_timecourse.samples, 1), 1);

NeuBlockFilt = strcmp(behav_data.Condition, 'neu');
NeuIDXs = block_durations_as_IDXs(NeuBlockFilt);

NegBlockFilt = strcmp(behav_data.Condition, 'neg');
NegIDXs = block_durations_as_IDXs(NegBlockFilt);

% Remove the overlapping elements
f = ismember(NeuIDXs, NegIDXs);
k = ismember(NegIDXs, NeuIDXs);
overlap = NeuIDXs(f);
NeuIDXs(f) = [];
NegIDXs(k) = [];

ds_timecourse.sa.EmotionalValence(NeuIDXs) = {'Neutral'};
ds_timecourse.sa.EmotionalValence(NegIDXs) = {'Negative'};

% ContextNumber = the context number that is occuring at each timepoint.
% n/a when fixation is on the screen.
ds_timecourse.sa.ContextNumber    = NaN(size(ds_timecourse.samples, 1), 1);

% force the ContextNum column of block onset events to take on the value 
% of the first trial within the miniblock
behav_data.ContextNum(isnan(behav_data.ContextNum)) = behav_data.ContextNum(find(isnan(behav_data.ContextNum)) + 1);

% A logical filter identifying block onset events
blockOnsetsFilt = isnan(behav_data.ColorNum);

Con1Filt = blockOnsetsFilt & behav_data.ContextNum == 1;
Con1IDXs = block_durations_as_IDXs(Con1Filt);

Con2Filt = blockOnsetsFilt & behav_data.ContextNum == 2;
Con2IDXs = block_durations_as_IDXs(Con2Filt);

Con3Filt = blockOnsetsFilt & behav_data.ContextNum == 3;
Con3IDXs = block_durations_as_IDXs(Con3Filt);

Con4Filt = blockOnsetsFilt & behav_data.ContextNum == 4;
Con4IDXs = block_durations_as_IDXs(Con4Filt);

Con1IDXs(ismember(Con1IDXs, overlap)) = [];
Con2IDXs(ismember(Con2IDXs, overlap)) = [];
Con3IDXs(ismember(Con3IDXs, overlap)) = [];
Con4IDXs(ismember(Con4IDXs, overlap)) = [];

ds_timecourse.sa.ContextNumber(Con1IDXs) = 1;
ds_timecourse.sa.ContextNumber(Con2IDXs) = 2;
ds_timecourse.sa.ContextNumber(Con3IDXs) = 3;
ds_timecourse.sa.ContextNumber(Con4IDXs) = 4;

ds_timecourse.sa.TimePoint = [1:size(ds_timecourse.samples, 1)]';

%% Calculate the informational timecourses for each ROI
% One informational timecourse per ROI. See cosmo_informational_timecourse

[ds_t, ds_a, ds_d] = cosmo_informational_timecourse(ds_template, ds_timecourse, 'EmotionalValence');

ds_t = cosmo_split(ds_t, {'ROIname'}, 2); % target info timecourse
ds_a = cosmo_split(ds_a, {'ROIname'}, 2); % alternative info timecourse
ds_d = cosmo_split(ds_d, {'ROIname'}, 2); % target - next highest alternative; "discrimination"

for ri = 1:numberOfROIs
        
    % Figure named after ROI
    figure('Name', char(ds_t{ri}.fa.ROIname), 'visible', 'off');
    
    % Target Condition Timecourse
    subplot(3, 1, 1)
    plot(ds_t{ri}.samples)
    title('Target MV Timecourse')
    axis([1 size(ds_t{ri}.samples, 1) -1 1])
    
    % Alternate Condition Timecourse
    subplot(3, 1, 2)
    plot(ds_a{ri}.samples)
    title('Alternate MV Timecourse(s)')
    axis([1 size(ds_a{ri}.samples, 1) -1 1])
    labels = cell(1, length(ds_a{1}.fa.AltCondNum));
    for f = 1:length(ds_a{1}.fa.AltCondNum)
        labels{f} = sprintf('Alternate Condition %d', ds_a{ri}.fa.AltCondNum(f));
    end
    legend(labels)
    
    % Discrimination Timecourse
    subplot(3, 1, 3)
    plot(ds_d{ri}.samples)
    title('Discrimination MV Timecourse')
    axis([1 size(ds_d{ri}.samples, 1) -1 1])
    
    % save
    filename = sprintf('%s_ROI-%s_information-%s_multivar-discim-timecourse.fig', subjects{s}, char(ds_t{ri}.fa.ROIname), information);
    savefig(fullfile(outpath, filename));
    
end

%% Correlate informational timecourses
% Calculate the correlation between each of the ROIs information
% timecoruses. Write out.

% the network, as a matrix
ds_d = cosmo_stack(ds_d, 2);
R = corrplot(ds_d.samples, 'type', 'Spearman', 'varNames', regexprep(ds_d.fa.ROIname, '_', ' '));
filename = sprintf('%s_network-%s_class-%s_corrplot.fig', subjects{s}, maskType, information);
savefig(fullfile(outpath, filename));

% the network, as a graph
figure;
G = graph(R, regexprep(ds_d.fa.ROIname, '_', ' '), 'OmitSelfLoops');
p = plot(G, 'LineWidth', abs(G.Edges.Weight)*100);
switch maskType
    case 'PM-System'
        % custom positions
        p.XData = [-1 1 -1 1];
        p.YData = [1 1 -1 -1];
    case 'AT-System'
        % custom positions
        p.XData = [-1 1 -1 1 -1.5 1.5];
        p.YData = [1 1 -1 -1 0 0];
end

% custom color code
colormap([1 0 0; 0 0 1])
CData = NaN(1, G.numedges);
CData(G.Edges.Weight < 0) = -1;
CData(G.Edges.Weight > 0) = 1;
p.EdgeCData = CData;
p.EdgeColor = 'flat';

str = sprintf('%s class-%s networkName-%s', subjects{s}, information, maskType);
title(str)
filename = sprintf('%s_network-%s_class-%s_graph.fig', subjects{s}, maskType, information);
savefig(fullfile(outpath, filename))

filename = sprintf('sub-%s_results.mat', subjects{s});
save(fullfile(outpath, filename), 'R', 'ds_t', 'ds_a', 'ds_d')

%% Subfunctions

function fullpathtofourDfile = threeDTofourD(files)

    fourD_filename = sprintf('sub-%s_res.nii', subjects{s});
    
    fullpathtofourDfile = fullfile(unique(cellfun(@fileparts, files, 'UniformOutput', false)), fourD_filename);
    
    if ~exist(char(fullpathtofourDfile), 'file')
    
        matlabbatch{1}.spm.util.cat.vols  = files;
        matlabbatch{1}.spm.util.cat.name  = fourD_filename;
        matlabbatch{1}.spm.util.cat.dtype = 0;

        spm_jobman('run', matlabbatch)
    
    end
    
end

function behav_data = custom_read_tsv(x)
    behav_data = readtable(x, 'FileType', 'text', 'Delimiter', '\t');
    behav_data.Run = repmat(regexp(x, 'run-0[0-6]', 'match'), height(behav_data), 1);
end

function IDXs = block_durations_as_IDXs(Filt)
    
    % Onsets and Offsets in TRs
    % Convert Onset and Offset times relative to beginning of run -->
    % relative to beginning of experiment   

    blockOnsetsInTRs  = [];
    blockOffsetsInTRs = [];
    
    runs = unique(behav_data.Run);
    for r = 1:length(runs)
        runFilt = strcmp(behav_data.Run, runs{r});
        F = runFilt & Filt;
        OffSet = (r - 1) * TRsPerSession;
        blockOnsetsInTRs  = vertcat(blockOnsetsInTRs, round(behav_data.onset(F, :) / TR + OffSet));
        blockOffsetsInTRs = vertcat(blockOffsetsInTRs, round((behav_data.onset(find(F) + 4) + 5) / TR + OffSet));
    end
    
    % Collect indexs
    IDXs = [];
    for b = 1:length(blockOnsetsInTRs)
        IDXs = horzcat(IDXs, blockOnsetsInTRs(b):blockOffsetsInTRs(b));
    end
end

end