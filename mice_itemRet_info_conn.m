function mice_itemRet_info_conn(ii)

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
outpath         = '/gsfs0/scratch/kurkela/results/mice-itemret-informational-connectivity';

%% Subjects
% Subject IDs
subjectsLevels           = cellstr(spm_select('List', data_path, 'dir', 's0[0-3][0-9]'));
analysesLevels           = {'wholeMTL'; 'withinMTL'};
class_information_Levels = {'ContextNum'}; % 'EmotionalValence';

%% Combinations of input parameters

% All combination of the input parmaters "subjects" and "analyses"
c=0;
for q = 1:length(class_information_Levels)
    for w = 1:length(analysesLevels)
        for r = 1:length(subjectsLevels)
            c=c+1;
            subject{c}     = subjectsLevels{r}; %#ok<*AGROW>
            analyses{c}    = analysesLevels{w};
            information{c} = class_information_Levels{q};
        end
    end
end

Combos = table(analyses', subject', information', ...
               'VariableNames', {'analysisType', 'subject', 'information'});

% Parameters for the current iteration
analysisType = Combos.analysisType{ii};
subject      = Combos.subject{ii};
information  = Combos.information{ii};

%% Masks
% full brain and MTL masks

switch analysisType
    
    case 'AT_PM_Systems'
        
        %%% PM-System
        
%             % left
%             masks{1}  = fullfile(rosies_roi_path, 'rHIPP_BODY_L_mask.nii');
%             masks{2}  = fullfile(rosies_roi_path, 'rPHC_ANT_L_mask.nii');
%             masks{3}  = fullfile(roses_roi_path, 'PCC_L_ROI.nii');
%             masks{4}  = fullfile(roses_roi_path, 'RSC_L_ROI.nii');
%             masks{5}  = fullfile(roses_roi_path, 'PREC_L_ROI.nii');
%             masks{6}  = fullfile(roses_roi_path, 'ANG_L_ROI.nii');
% 
%             % right
%             masks{7}  = fullfile(rosies_roi_path, 'rHIPP_BODY_R_mask.nii');
%             masks{8}  = fullfile(rosies_roi_path, 'rPHC_ANT_R_mask.nii');
%             masks{9}  = fullfile(roses_roi_path, 'PCC_R_ROI.nii');
%             masks{10} = fullfile(roses_roi_path, 'RSC_R_ROI.nii');
%             masks{11} = fullfile(roses_roi_path, 'PREC_R_ROI.nii');
%             masks{12} = fullfile(roses_roi_path, 'ANG_R_ROI.nii');
% 
        %%% AT-System
% 
%             % left
%             masks{13} = fullfile(rosies_roi_path, 'rAMY_L_mask.nii');
%             masks{14} = fullfile(rosies_roi_path, 'rHIPP_HEAD_L_mask.nii');
%             masks{15} = fullfile(rosies_roi_path, 'rPRC_L_mask.nii');
%             masks{16} = fullfile(roses_roi_path, 'ITC_L_ROI.nii');
%             masks{17} = fullfile(roses_roi_path, 'OFC_L_ROI.nii');
% 
%             % right
%             masks{18} = fullfile(rosies_roi_path, 'rAMY_R_mask.nii');
%             masks{19} = fullfile(rosies_roi_path, 'rHIPP_HEAD_R_mask.nii');
%             masks{20} = fullfile(rosies_roi_path, 'rPRC_R_mask.nii');
%             masks{21} = fullfile(roses_roi_path, 'ITC_R_ROI.nii');
%             masks{22} = fullfile(roses_roi_path, 'OFC_R_ROI.nii');
        
       ROIs   = {'HIPP_BODY_L', ...
                 'PHC_ANT_L', ...
                 'PCC_L_ROI', ...
                 'PCC_L_ROI', ...
                 'PREC_L_ROI', ...
                 'ANG_L_ROI', ...
                 'HIPP_BODY_R', ...
                 'PHC_ANT_R', ...
                 'PCC_R_ROI', ...
                 'RSC_R_ROI', ...
                 'PREC_R_ROI', ...
                 'ANG_R_ROI', ...
                 'AMY_L', ...
                 'HIPP_HEAD_L', ...
                 'PRC_L', ...
                 'ITC_L_ROI', ...
                 'OFC_L_ROI', ...
                 'AMY_R', ...
                 'HIPP_HEAD_R', ...
                 'PRC_R', ...
                 'ITC_R_ROI', ...
                 'OFC_R_ROI'};
                     
    case 'wholeMTL'

        %%% MTL

%             masks{1}  = fullfile(outpath, 'rois', 'rMTL_group50_mask_L.nii');
%             masks{2}  = fullfile(outpath, 'rois', 'rMTL_group50_mask_R.nii');
% 
%         %%% PM-System
% 
%             % left
%             masks{3}  = fullfile(roses_roi_path, 'PCC_L_ROI.nii');
%             masks{4}  = fullfile(roses_roi_path, 'RSC_L_ROI.nii');
%             masks{5}  = fullfile(roses_roi_path, 'PREC_L_ROI.nii');
%             masks{6}  = fullfile(roses_roi_path, 'ANG_L_ROI.nii');
% 
%             % right
%             masks{7}  = fullfile(roses_roi_path, 'PCC_R_ROI.nii');
%             masks{8} = fullfile(roses_roi_path, 'RSC_R_ROI.nii');
%             masks{9} = fullfile(roses_roi_path, 'PREC_R_ROI.nii');
%             masks{10} = fullfile(roses_roi_path, 'ANG_R_ROI.nii');
% 
%         %%% AT-System
% 
%             % left
%             masks{11} = fullfile(roses_roi_path, 'ITC_L_ROI.nii');
%             masks{12} = fullfile(roses_roi_path, 'OFC_L_ROI.nii');
% 
%             % right
%             masks{13} = fullfile(roses_roi_path, 'ITC_R_ROI.nii');
%             masks{14} = fullfile(roses_roi_path, 'OFC_R_ROI.nii');

        ROIs = {'MTL_group50_L', ...
                'MTL_group50_R', ...
                'PCC_L_ROI', ...
                'RSC_L_ROI', ...
                'PREC_L_ROI', ...
                'ANG_L_ROI', ...
                'PCC_R_ROI', ...
                'RSC_R_ROI', ...
                'PREC_R_ROI', ...
                'ANG_R_ROI', ...
                'ITC_L_ROI', ...
                'OFC_L_ROI', ...
                'ITC_R_ROI', ...
                'OFC_R_ROI'};

    case 'withinMTL'
        
        %%% PM-System
        
%             % left
%             masks{1}  = fullfile(rosies_roi_path, 'rHIPP_BODY_L_mask.nii');
%             masks{2}  = fullfile(rosies_roi_path, 'rPHC_ANT_L_mask.nii');
%             
%             % right
%             masks{3}  = fullfile(rosies_roi_path, 'rHIPP_BODY_R_mask.nii');
%             masks{4}  = fullfile(rosies_roi_path, 'rPHC_ANT_R_mask.nii');
%             
%         %%% AT-System
%             
%             % left
%             masks{5} = fullfile(rosies_roi_path, 'rAMY_L_mask.nii');
%             masks{6} = fullfile(rosies_roi_path, 'rHIPP_HEAD_L_mask.nii');
%             masks{7} = fullfile(rosies_roi_path, 'rPRC_L_mask.nii');
%             
%             % right
%             masks{8} = fullfile(rosies_roi_path, 'rAMY_R_mask.nii');
%             masks{9} = fullfile(rosies_roi_path, 'rHIPP_HEAD_R_mask.nii');
%             masks{10} = fullfile(rosies_roi_path, 'rPRC_R_mask.nii');
            
        ROIs = {'HIPP_BODY_L', ...
                'PHC_ANT_L', ...
                'HIPP_BODY_R', ...
                'PHC_ANT_R', ...
                'AMY_L', ...
                'HIPP_HEAD_L', ...
                'PRC_L', ...
                'AMY_R', ...
                'HIPP_HEAD_R', ...
                'PRC_R'};
        
end

%% Load Single Trial Estimates

ds = [];
load(fullfile(outpath, 'SingleTrialEstimates', sprintf('sub-%s_SingleTrialEstimates.mat', subject)))

ds = cosmo_slice(ds, ismember(ds.fa.ROIlabel, ROIs), 2);   
ds = cosmo_slice(ds, ~strcmp(ds.sa.Probe, 'novel'));

outpath = fullfile(outpath, analysisType);
if ~exist(outpath, 'dir')
    mkdir(outpath)
end

%% Parameters
% hemodynamic_lag (in TRs).
numberOfROIs        = length(unique(ds.fa.ROIlabel));

%% Timecourses
% Instead of a "timecourse", we are going to correlates across a beta
% series. To hack this setup, just set ds_timecourse = ds_template and add
% a "TimePoint" sample attribute.

ds_timecourse = ds;
ds_timecourse.sa.TimePoint = (1:size(ds_timecourse.samples, 1))';

%% Calculate the informational timecourses for each ROI
% One informational timecourse per ROI. See cosmo_informational_timecourse

%%% Actual FC

[ds_t, ds_a, ds_d] = cosmo_informational_timecourse(ds, ds_timecourse, information);

%%% Null FC

%-- Find balanced random permutations
% The startegy here is brute force: randomly permutate the .(information)
% field until we find a permutation that balances the distribution factor
% levels. Note: if there are many factor levels, this code if HIGHLY
% INEFFICIENT.
nSims        = 1000;
uniqueChunks = unique(ds.sa.chunks);
tmp          = tabulate(ds.sa.chunks);
GoodPerms    = nan(unique(tmp(:, 2)), length(uniqueChunks), nSims);

for c = 1:length(uniqueChunks)

    % The current chunk's indices, as a vector; the target trial labels as
    % they are currently organized; initalizing count variables for the
    % while loop.
    thisChunksIdxs         = find(ds.sa.chunks == uniqueChunks(c));
    oldTargets             = ds.sa.(information)(thisChunksIdxs);
    j = 0; i = 0;
    
    while i < nSims
        
        % The iteration number, a random permutation of this chunks's
        % indices, the oldTargets, the newTargets, and a crosstabulation
        % of oldTargets x newTargets
        j = j+1;
        thisChunksIdxsRandPerm = thisChunksIdxs(randperm(length(thisChunksIdxs)));
        newTargets             = ds.sa.(information)(thisChunksIdxsRandPerm);
        tbl                    = crosstab(oldTargets, newTargets);
        
        % Balanced Tabulate Table
        % Create a matrix that will force the tabluated random permuation
        % to be--roughly speaking--balanced.
        Counts               = tabulate(oldTargets);
        numberOfFactorLevels = size(Counts, 1);
        if isa(Counts, 'cell')
            balanceTabulate = repmat(cell2mat(Counts(:, 2)) / numberOfFactorLevels, 1, numberOfFactorLevels);
        elseif isa(Counts, 'double')
            balanceTabulate = repmat(Counts(:, 2) / numberOfFactorLevels, 1, numberOfFactorLevels);
        end
        
        % Gather up the "Good" permutations that are balanced across the
        % factor levels.
        if all(ismembertol(tbl, balanceTabulate, 1, 'DataScale', 1))
           i = i + 1;
           fprintf('Good Perm %d ... \n\n', i);
           GoodPerms(1:length(thisChunksIdxsRandPerm), c, i) = thisChunksIdxsRandPerm;
        end
        
    end
end

% Run Null Simulations
for e = 1:size(GoodPerms, 3)

    fprintf('Simulation %d...\n\n', e)
    
    % a copy of the original data set
    ds_null = ds;
    
    % a vector of randomly permutated indices that are balanced across
    % levels of the factor "information"
    lengthOfTheThisSimsGoodPermIDxsVector = size(GoodPerms, 1) * size(GoodPerms, 2);
    thisSimsGoodPermIDxs = reshape(GoodPerms(:,:,e), lengthOfTheThisSimsGoodPermIDxsVector, 1);
    
    % Randomly relabeling the .information field in a balanced manner
    ds_null.sa.(information) = ds_null.sa.(information)(thisSimsGoodPermIDxs);
    
    % Run a null simulation!
    [~, ~, null_ds_d] = cosmo_informational_timecourse(ds_null, ds_timecourse, information);

    % The Null Correlation Matrix
    if e == 1
        nullR = corr(null_ds_d.samples);
    else
        nullR = cat(3, nullR, corr(null_ds_d.samples));
    end
    
end

figure;imagesc(mean(nullR, 3)); 
cb1 = colorbar; cb1.Label.String = 'mean r';
title(sprintf('Mean Null Corr Matrix After %d Simulations', e))
[minC(1), maxC(1)] = caxis();
axis1 = gca;

figure;imagesc(corr(ds_d.samples)); 
cb2 = colorbar; cb2.Label.String = 'r';
title('Corr Matrix')
[minC(2), maxC(2)] = caxis();
axis2 = gca;

figure;imagesc(corr(ds_d.samples) - mean(nullR, 3)); 
cb3 = colorbar; cb3.Label.String = 'Difference';
title('Corr Matrix - Null Corr Matrix')
[minC(3), maxC(3)] = caxis();
axis3 = gca;
caxis(axis1, [min(minC) max(maxC)])
caxis(axis2, [min(minC) max(maxC)])


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
%     filename = sprintf('%s_ROI-%s_information-%s_multivar-discim-timecourse.fig', subject, char(ds_t{ri}.fa.ROIname), information);
%     savefig(fullfile(outpath, filename));
    
end

%% Correlate informational timecourses
% Calculate the correlation between each of the ROIs information
% timecoruses. Write out.

% the network, as a matrix
ds_d = cosmo_stack(ds_d, 2);
R = corrplot(ds_d.samples, 'type', 'Spearman', 'varNames', regexprep(ds_d.fa.ROIname, '_', ' '));
% filename = sprintf('%s_class-%s_corrplot.fig', subject, information);
% savefig(fullfile(outpath, filename));

% the network, as a graph
figure;
G = graph(R, regexprep(ds_d.fa.ROIname, '_', ' '), 'OmitSelfLoops');
p = plot(G, 'LineWidth', abs(G.Edges.Weight)*100);

% custom color code
colormap([1 0 0; 0 0 1])
CData = NaN(1, G.numedges);
CData(G.Edges.Weight < 0) = -1;
CData(G.Edges.Weight > 0) = 1;
p.EdgeCData = CData;
p.EdgeColor = 'flat';

str = sprintf('%s class-%s', subject, information);
title(str)
% filename = sprintf('%s_class-%s_graph.fig', subject, information);
% savefig(fullfile(outpath, filename))

filename = sprintf('sub-%s_class-%s_results.mat', subject, information);
save(fullfile(outpath, filename), 'R', 'nullR', 'ds_t', 'ds_a', 'ds_d', 'ds')

end
