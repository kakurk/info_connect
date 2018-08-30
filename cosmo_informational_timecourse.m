function [ discriminability_MV_timecourse ] = cosmo_informational_timecourse( ds_betas, ds_timecourse, targetCondition, value )
% Using a leave-one-run-out procedure, return an "informational timecourse"
%

%% Parse input arguments

assert(isequal(unique(ds_betas.sa.chunks), unique(ds_timecourse.sa.chunks)), 'both datasets must have the same number of chunks')

chunks     = unique(ds_betas.sa.chunks);
discriminability_MV_timecourse = cell(1, length(chunks));

for ch = chunks'
    
    % Template Patterns
    Filt = ~ismember(ds_betas.sa.chunks, ch);
    ds_template_patterns = cosmo_slice(ds_betas, Filt);
    ds_template_patterns = cosmo_average_samples(ds_template_patterns, 'split_by', {targetCondition});
    
    % Timecourse
    Filt = ismember(ds_timecourse.sa.chunks, ch);
    ds_multivariate_timecourse = cosmo_slice(ds_timecourse, Filt);
    
    % Target Pattern MV Timecourse
    if iscellstr(ds_template_patterns.sa.(targetCondition))
        Filt              = strcmp(ds_template_patterns.sa.(targetCondition), value);
    elseif isnumeric(ds_template_patterns.sa.(targetCondition))
        Filt              = ds_template_patterns.sa.(targetCondition) == value;
    end
    ds_target_pattern = cosmo_slice(ds_template_patterns, Filt);
    target_pattern_MV_timecourse = cosmo_corr(ds_target_pattern.samples', ds_multivariate_timecourse.samples');
    target_pattern_MV_timecourse = atanh(target_pattern_MV_timecourse);
    
    % Next Highest Alternate Pattern
    if iscellstr(ds_template_patterns.sa.(targetCondition))
        Filt              = ~strcmp(ds_template_patterns.sa.(targetCondition), value);
    elseif isnumeric(ds_template_patterns.sa.(targetCondition))
        Filt              = ds_template_patterns.sa.(targetCondition) ~= value;
    end
    ds_target_pattern = cosmo_slice(ds_template_patterns, Filt);
    alternate_pattern_MV_timecourse = cosmo_corr(ds_target_pattern.samples', ds_multivariate_timecourse.samples');
    alternate_pattern_MV_timecourse = max(alternate_pattern_MV_timecourse, [], 1);
    alternate_pattern_MV_timecourse = atanh(alternate_pattern_MV_timecourse);
    
    % 'Discriminability'
    discriminability_MV_timecourse{ch} = target_pattern_MV_timecourse - alternate_pattern_MV_timecourse;
    
end
    
discriminability_MV_timecourse = horzcat(discriminability_MV_timecourse{:});

end