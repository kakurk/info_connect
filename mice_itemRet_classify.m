function mice_itemRet_classify(s)
% mice_itemRet_classify. Calculate mean classification accuracy for each 
% node.

%% Add Paths
% SPM12, CoSMoMVPA, and the Informational Connectivity Toolbox
addpath(genpath('/gsfs0/data/kurkela/Documents/toolboxes-fmri/spm12'));
addpath(genpath('/gsfs0/data/kurkela/Documents/toolboxes-fmri/CoSMoMVPA'));
addpath(genpath('/gsfs0/data/kurkela/Documents/libsvm'));

%% Relevant Directories
% data_path = where the single trial beta images are
% inital_spm_model_path = where the original SPM models are
% bids_path = where the raw BIDS formatted data are
data_path       = '/gsfs0/data/ritcheym/data/fmri/mice/analysis/retrieval/model_est/lssGenBeta-output/unsmoothed';
outpath         = '/gsfs0/scratch/kurkela/results/mice-itemret-informational-connectivity';

outpath = fullfile(outpath, 'classification_at_nodes');
if ~exist(outpath, 'dir')
    mkdir(outpath)
end

%% Subjects
% Subject IDs
subjects = cellstr(spm_select('List', data_path, 'dir', 's0[0-3][0-9]'));

%% Load Single Trial Estimates
ds = [];
load(fullfile(fileparts(outpath), 'SingleTrialEstimates', sprintf('sub-%s_SingleTrialEstimates.mat', subjects{s})))

%% Classify!

% initalize a targets sample attribute. establish what we are classifying
ds.sa.targets = NaN(size(ds.samples, 1), 1);
what_we_are_classifying = 'ContextNum';

% define the targets sample attribute based on what we want to classify
switch what_we_are_classifying
    case 'EmotionalValence'
        NegativeTrialsFilter = strcmp(ds.sa.EmotionalValence, 'neg');
        ds.sa.targets(NegativeTrialsFilter) = 1;
        NeutralTrialsFilter = strcmp(ds.sa.EmotionalValence, 'neut');
        ds.sa.targets(NeutralTrialsFilter) = 2;
    case 'ContextNum'
        ds.sa.targets(ds.sa.ContextNum == 1) = 1;
        ds.sa.targets(ds.sa.ContextNum == 2) = 2;
        ds.sa.targets(ds.sa.ContextNum == 3) = 3;
        ds.sa.targets(ds.sa.ContextNum == 4) = 4;
end

% Slice out any NaN targets
ds = cosmo_slice(ds, ~isnan(ds.sa.targets));

% partition
partitions = cosmo_nfold_partitioner(ds);
partitions = cosmo_balance_partitions(partitions, ds);
fprintf('There are %d partitions\n', numel(partitions.train_indices));
fprintf('# train samples:%s\n', sprintf(' %d', cellfun(@numel, ...
    partitions.train_indices)));
fprintf('# test samples:%s\n', sprintf(' %d', cellfun(@numel, ...
    partitions.test_indices)));

ds_roi = cosmo_split(ds, {'ROIlabel'}, 2);

for i = 1:length(ds_roi)

    % Run Classification
    [~, accuracy] = cosmo_crossvalidate(ds_roi{i}, ...
                                        @cosmo_classify_libsvm, ...
                                        partitions);

    % Report Results
    fprintf('Accuracy: %0.2f\n', accuracy);
    
    % Record Results
    subject = cellstr(subjects{s});
    roi     = cellstr(unique(ds_roi{i}.fa.ROIlabel));
    class   = cellstr(what_we_are_classifying);
    if i == 1
        results = table(subject, roi, class, accuracy);
    else
        results = vertcat(results, table(subject, roi, class, accuracy));
    end
    
end

% Define output location
output_fn = fullfile(outpath, sprintf('sub-%s_class-%s_svmClassificationResults.csv', subjects{s}, what_we_are_classifying));

% Store results to disc
writetable(results, output_fn)

end