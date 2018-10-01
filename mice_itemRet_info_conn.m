function mice_itemRet_info_conn(s)

%% Add Paths
% SPM12, CoSMoMVPA, and the Informational Connectivity Toolbox
addpath(genpath('/gsfs0/data/kurkela/Documents/toolboxes-fmri/spm12'));
addpath(genpath('/gsfs0/data/kurkela/Documents/toolboxes-fmri/CoSMoMVPA'));
addpath(genpath('/gsfs0/data/kurkela/Documents/info_conn/thirdparty'));

%% Relevant Directories
% data_path = where the single trial beta images are
% inital_spm_model_path = where the original SPM models are
% bids_path = where the raw BIDS formatted data are
data_path  = '/gsfs0/data/ritcheym/data/fmri/mice/analysis/retrieval/model_est/lssGenBeta-output/unsmoothed';
roi_path   = '/gsfs0/data/ritcheym/data/fmri/mice/analysis/retrieval/model_est/lssGenBeta-output/unsmoothed/ROIs';
outpath    = '/gsfs0/scratch/kurkela/results/mice-itemret-informational-connectivity';

%% Masks
% full brain and MTL masks
maskType = 'PM-System';
switch maskType
    case 'PM-System'
        % left
        masks{1} = fullfile(roi_path, 'rHIPP_BODY_L_mask.nii');
        masks{2} = fullfile(roi_path, 'rPHC_ANT_L_mask.nii');
        % right
        masks{3} = fullfile(roi_path, 'rHIPP_BODY_R_mask.nii'); 
        masks{4} = fullfile(roi_path, 'rPHC_ANT_R_mask.nii');
    case 'AT-System'
        % left
        masks{1} = fullfile(roi_path, 'rAMY_L_mask.nii');
        masks{2} = fullfile(roi_path, 'rHIPP_HEAD_L_mask.nii');
        masks{3} = fullfile(roi_path, 'rPRC_L_mask.nii');
        % right
        masks{4} = fullfile(roi_path, 'rAMY_R_mask.nii');
        masks{5} = fullfile(roi_path, 'rHIPP_HEAD_R_mask.nii');
        masks{6} = fullfile(roi_path, 'rPRC_R_mask.nii');
end

%% Subjects
% Subject IDs
subjects = cellstr(spm_select('List', data_path, 'dir', 's0[0-3][0-9]'));

%% Parameters
% hemodynamic_lag (in TRs).
numberOfROIs        = length(masks);
information         = 'EmotionalValence'; % 'ContextNumber'

%% Templates
% In order to perform informational connectivity, we need to define several
% "template patterns" to correlate with the fMRI timecourse. There should
% be one "template pattern" per condition, that is defined based on
% independent training data.

% all of the single trial beta images for this subject
single_trial_betasFNs = RecurseAndFilterFileSearch(data_path, 'Sess00[1-6].*\.nii', [filesep subjects{s}]);
single_trial_betasFNs(~contains(single_trial_betasFNs, 'itemret')) = [];

% We need to 4-d file them.
fourDsingletrialbetasFN = threeDTofourD(single_trial_betasFNs, outpath);

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
regularExpression               = '(?<=_)ne[gu][t]?(?=-)';
matches                         = regexp(single_trial_betasFNs, regularExpression, 'match');
matches                         = cellfun(@char, matches, 'UniformOutput', false);
ds_template.sa.EmotionalValence = matches;

% using the custom regularExpression, extract the substring from the
% single_trial_betaFNs identifying the beta's emotional valence. Take
% the match and convert it from a nested cell array to a cell string.
% Assign it as a new sample attribute "EmotionalValence".
regularExpression    = '[tln][auo][rv][ge][l]?(?=-)';
matches              = regexp(single_trial_betasFNs, regularExpression, 'match');
matches              = cellfun(@char, matches, 'UniformOutput', false);
ds_template.sa.Probe = matches;

% using the custom regularExpression, extract the substring from the
% single_trial_betaFNs identifying the beta's emotional valence. Take
% the match and convert it from a nested cell array to a cell string.
% Assign it as a new sample attribute "EmotionalValence".
regularExpression     = '(?<=-)[hmCF][iRA][ts]?[s]?(?=_)';
matches               = regexp(single_trial_betasFNs, regularExpression, 'match');
matches               = cellfun(@char, matches, 'UniformOutput', false);
ds_template.sa.Memory = matches;

% using custom regularExpression, extract the session number from the
% single_trial_betaFNs.
regularExpression     = '(?<=Sess00)[0-6](?=_)';
matches               = regexp(single_trial_betasFNs, regularExpression, 'match');
matches               = cellfun(@(x) str2double(cell2mat(x)), matches);
ds_template.sa.chunks = matches;

% And, finally, the context information...which is not contained within the
% beta file names :(.
singleTrialInfoMat = RecurseAndFilterFileSearch(data_path, 'MICE_singletrial.*\.mat', subjects{s});
load(char(singleTrialInfoMat), 'subTable');
function out = custom_extraction(x)
    if ~isempty(x)
        out = cell2mat(x);
    else
        out = NaN;
    end
end
ds_template.sa.ContextNum = cellfun(@custom_extraction, subTable.context(contains(subTable.trial, 'itemret')));

%% Timecourses
% Instead of a "timecourse", we are going to correlates across a beta
% series. To hack this setup, just set ds_timecourse = ds_template and add
% a "TimePoint" sample attribute.

ds_timecourse = ds_template;
ds_timecourse.sa.TimePoint = (1:size(ds_timecourse.samples, 1))';

%% Calculate the informational timecourses for each ROI
% One informational timecourse per ROI. See cosmo_informational_timecourse

[ds_t, ds_a, ds_d] = cosmo_informational_timecourse(ds_template, ds_timecourse, 'EmotionalValence');

ds_t = cosmo_split(ds_t, {'ROIname'}, 2); % target info timecourse
ds_a = cosmo_split(ds_a, {'ROIname'}, 2); % alternative info timecourse
ds_d = cosmo_split(ds_d, {'ROIname'}, 2); % target - next highest alternative; "discrimination"

for ri = 1:numberOfROIs
        
    % Figure named after ROI
    figure('Name', char(ds_t{ri}.fa.ROIname), 'visible', 'on');
    
    % Target Condition Timecourse
    subplot(3, 1, 1)
    plot(ds_t{ri}.samples)
    title('Target MV Timecourse')
    xticks(1:32:size(ds_t{ri}.samples, 1))
    axis([1 size(ds_t{ri}.samples, 1) -1 1])
    
    % Alternate Condition Timecourse
    subplot(3, 1, 2)
    plot(ds_a{ri}.samples, 'g')
    title('Alternate MV Timecourse(s)')
    axis([1 size(ds_a{ri}.samples, 1) -1 1])
    xticks(1:32:size(ds_a{ri}.samples, 1))
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
    xticks(1:32:size(ds_d{ri}.samples, 1))
    
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