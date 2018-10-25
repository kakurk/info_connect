function mice_itemRet_univariate(s)
% mice_itemRet_univariate. Calculate mean activation within each ROI

%% Add Paths
% SPM12, CoSMoMVPA, and the Informational Connectivity Toolbox
addpath(genpath('/gsfs0/data/kurkela/Documents/toolboxes-fmri/spm12'));
addpath(genpath('/gsfs0/data/kurkela/Documents/toolboxes-fmri/CoSMoMVPA'));

%% Relevant Directories
% data_path = where the single trial beta images are
% inital_spm_model_path = where the original SPM models are
% bids_path = where the raw BIDS formatted data are
data_path       = '/gsfs0/data/ritcheym/data/fmri/mice/analysis/retrieval/model_est/lssGenBeta-output/unsmoothed';
outpath         = '/gsfs0/scratch/kurkela/results/mice-itemret-informational-connectivity/';

outpath = fullfile(outpath, 'univariate_at_nodes');
if ~exist(outpath, 'dir')
    mkdir(outpath)
end

%% Subjects
% Subject IDs
subjects = cellstr(spm_select('List', data_path, 'dir', 's0[0-3][0-9]'));

%% Load Single Trial Estimates

ds = [];
load(fullfile(fileparts(outpath), 'SingleTrialEstimates', sprintf('sub-%s_SingleTrialEstimates.mat', subjects{s})))

%% Calculate Univariate Activation

% take the mean across voxels within each ROI
ds_mean   = cosmo_fx(ds, @(x) mean(x, 2), {'ROIlabel'}, 2);

% convert the cosmo data structure --> tidyverse long table. see the
% subfunction below for the specifics
results = cosmo2tidyverse(ds_mean);

% remove
results.i = [];
results.j = [];
results.k = [];

% rename
results.Properties.VariableNames{contains(results.Properties.VariableNames, 'samples')} = 'meanActivation';

% resort
results = results(:, {'ExperimentName', 'SubjectID', 'chunks', 'ContextNum', 'EmotionalValence', 'Memory', 'Probe', 'ROIlabel', 'meanActivation'});

% Define output location
output_fn = fullfile(outpath, sprintf('sub-%s_univariateResults.csv', subjects{s}));

% Store results to disc
writetable(results, output_fn)

%% Subfunctions

function tidyverseTable = cosmo2tidyverse(cosmo_ds)

    % Extract the samples into one long vector
    samples    = num2cell(cosmo_ds.samples(:));

    % match the feature attributes
    featureAtt = fieldnames(cosmo_ds.fa);
    FAtableCells = cell(length(samples), length(featureAtt));
    for f = 1:length(featureAtt)
        if isa(cosmo_ds.fa.(featureAtt{f}), 'double')
            FAtableCells(:,f) = num2cell(repelem(cosmo_ds.fa.(featureAtt{f})', size(cosmo_ds.samples,1)));
        elseif isa(cosmo_ds.fa.(featureAtt{f}), 'cell')
            FAtableCells(:,f) = repelem(cosmo_ds.fa.(featureAtt{f})', size(cosmo_ds.samples,1));
        end
    end

    % match the sample attributes
    sampleAtt    = fieldnames(cosmo_ds.sa);
    SAtableCells = cell(length(samples), length(sampleAtt));
    for f = 1:length(sampleAtt)
        if isa(cosmo_ds.sa.(sampleAtt{f}), 'double')
            SAtableCells(:,f) = num2cell(repmat(cosmo_ds.sa.(sampleAtt{f}), size(cosmo_ds.samples, 2), 1));
        elseif isa(cosmo_ds.sa.(sampleAtt{f}), 'cell')
            SAtableCells(:,f) = repmat(cosmo_ds.sa.(sampleAtt{f}), size(cosmo_ds.samples, 2), 1);
        end
    end

    % construct the tidyverse table
    tidyverseTable = horzcat(SAtableCells, FAtableCells, samples);
    tidyverseTable = cell2table(tidyverseTable, 'VariableNames', horzcat(sampleAtt', featureAtt', {'samples'}));

    % add overall attributes
    attributes = fieldnames(cosmo_ds.a);
    for a = 1:length(attributes)
        if ~isa(cosmo_ds.a.(attributes{a}), 'struct')
            tidyverseTable.(attributes{a}) = repmat(cosmo_ds.a.(attributes{a}), length(samples), 1);
        end
    end
    
end

end