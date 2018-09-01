function [ ds_t, ds_a, ds_d ] = cosmo_informational_timecourse( ds_betas, ds_timecourse, targetCondition )
% Using a leave-one-run-out procedure, return an "informational timecourse"
%

%% Parse input arguments

assert(isequal(unique(ds_betas.sa.chunks), unique(ds_timecourse.sa.chunks)), 'both datasets must have the same number of chunks')

if iscellstr(ds_timecourse.sa.(targetCondition))
    Filt = ~strcmp(ds_timecourse.sa.(targetCondition), 'n/a');
elseif isnumeric(ds_timecourse.sa.(targetCondition))
    Filt = ~isnan(ds_timecourse.sa.(targetCondition));
end
ds_timecourse = cosmo_slice(ds_timecourse, Filt);

chunks     = unique(ds_betas.sa.chunks);
ROIs       = unique(ds_betas.fa.ROIlabel);

for r = 1:length(ROIs)

    % current ROI name and filter
    currentROI = ROIs{r};
    currentROIFilt = strcmp(ds_betas.fa.ROIlabel, currentROI);
    
    % Slice off just the data from this ROI
    this_roi_ds_betas = cosmo_slice(ds_betas, currentROIFilt, 2);
    this_roi_ds_timecourse = cosmo_slice(ds_timecourse, currentROIFilt, 2);
    
    % the different factor levels
    condition_values = unique(ds_betas.sa.(targetCondition));
    
    % i.e., number of combinations of conditions and chunks. initalize the
    % proper number of cell arrarys
    ncombos        = length(condition_values) * max(chunks);
    target         = cell(1, ncombos);
    alternatives   = cell(1, ncombos);
    count          = 0;
    
    
    for ch = chunks'

        % Template Patterns NOT FROM THIS CHUNK
        Filt = ~ismember(this_roi_ds_betas.sa.chunks, ch);
        ds_template_patterns = cosmo_slice(this_roi_ds_betas, Filt);
        ds_template_patterns = cosmo_average_samples(ds_template_patterns, 'split_by', {targetCondition});

        % Timecourse FROM THIS CHUNK
        Filt = ismember(this_roi_ds_timecourse.sa.chunks, ch);
        ds_multivariate_timecourse = cosmo_slice(this_roi_ds_timecourse, Filt);
        
        for c = 1:length(condition_values)
            
            % advance the counter
            count = count + 1;
            
            % a string; current condition name
            currentCondition = condition_values(c);
            if iscellstr(currentCondition)
                currentCondition = char(currentCondition);
            end
            
            % The Target Timecourse
            if iscellstr(ds_multivariate_timecourse.sa.(targetCondition))
                Filt              = strcmp(ds_multivariate_timecourse.sa.(targetCondition), currentCondition);
            elseif isnumeric(ds_multivariate_timecourse.sa.(targetCondition))
                Filt              = ds_multivariate_timecourse.sa.(targetCondition) == currentCondition;
            end
            ds_multivariate_timecourse_c = cosmo_slice(ds_multivariate_timecourse, Filt);
        
            % The Target Pattern
            if iscellstr(ds_template_patterns.sa.(targetCondition))
                Filt              = strcmp(ds_template_patterns.sa.(targetCondition), currentCondition);
            elseif isnumeric(ds_template_patterns.sa.(targetCondition))
                Filt              = ds_template_patterns.sa.(targetCondition) == currentCondition;
            end
            ds_target_pattern_c = cosmo_slice(ds_template_patterns, Filt);
                   
            target{count}.sa      = ds_multivariate_timecourse_c.sa;
            target{count}.samples = cosmo_corr(ds_target_pattern_c.samples', ds_multivariate_timecourse_c.samples')';
            target{count}.samples = atanh(target{count}.samples);

            % The Alternate Patterns
            if iscellstr(ds_template_patterns.sa.(targetCondition))
                Filt              = ~strcmp(ds_template_patterns.sa.(targetCondition), currentCondition);
            elseif isnumeric(ds_template_patterns.sa.(targetCondition))
                Filt              = ds_template_patterns.sa.(targetCondition) ~= currentCondition;
            end
            ds_target_pattern_c = cosmo_slice(ds_template_patterns, Filt);
                        
            alternatives{count}.sa      = ds_multivariate_timecourse_c.sa;
            alternatives{count}.samples = cosmo_corr(ds_target_pattern_c.samples', ds_multivariate_timecourse_c.samples')';
            alternatives{count}.samples = atanh(alternatives{count}.samples);

        end
        
    end
    
    % stack
    target            = cosmo_stack(target);
    target.fa.ROIname = {currentROI};
    % reorder
    [~, sortIDxs]  = sort(target.sa.TimePoint);
    target.samples = target.samples(sortIDxs);
    target.sa      = structfun(@(x) x(sortIDxs), target.sa, 'UniformOutput', false);

    % stack
    alternatives               = cosmo_stack(alternatives);
    alternatives.fa.ROIname    = repmat({currentROI}, 1, size(alternatives.samples, 2));
    alternatives.fa.AltCondNum = 1:size(alternatives.samples, 2);
    
    % reorder
    [~, sortIDXs]        = sort(alternatives.sa.TimePoint);
    alternatives.samples = alternatives.samples(sortIDXs, :);
    alternatives.sa      = structfun(@(x) x(sortIDXs), alternatives.sa, 'UniformOutput', false);

    % discimination
    discrimination         = target; % initalize
    discrimination.samples = target.samples - max(alternatives.samples, [], 2);
    
    if r == 1
        ds_t = target;
        ds_a = alternatives;
        ds_d = discrimination;
    else
        ds_t = cosmo_stack({ds_t target}, 2);
        ds_a = cosmo_stack({ds_a alternatives}, 2);
        ds_d = cosmo_stack({ds_d discrimination}, 2);
    end
    
end

end