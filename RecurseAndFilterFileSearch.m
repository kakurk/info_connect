function files = RecurseAndFilterFileSearch(directory, regularExpression, Filter)
% Recursively search a directory for files and then filter the matches
%
% directory = fullpath string to the directory you would like to
% recursively search. Ex: '/users/kylekurkela/dataset/MICE'
%
% regularExpression = string identifying a regular expression used to
% search the directory. Ex: '.*\.nii'. AKA all .nii files.
%
% filter = string identifying a substring used to filter the results by.
% Ex: 'sub-s002'. Only grab files that have 'sub-s002' somewhere in the
% full path. This can be in either the filename or one of the
% subdirectories.
%
% files  =  cellstring containing the full paths to identified files. Ex:
% {'/users/kylekurkela/dataset/MICE/sub-s002/spmT_0001.nii',
%  '/users/kylekurkela/dataset/MICE/sub-s002/con_0001.nii'}
%% Parse Input Arguments

% directory
assert(ischar(directory), 'directory must be a string')
assert(exist(directory, 'dir') == 7, 'directory does not exist')

% regularExpression
assert(ischar(regularExpression), 'regularExpression must be a string')

% filter
assert(ischar(Filter), 'Filter must be a string')

%% Check for dependency

assert(~isempty(which('spm_select')), 'SPM must be on the matlab searchpath')

%% Recurse
% Using SPM's spm_select function, recursively select all files in the
% directory that match the regular expression. Throw an informative error
% if we don't find anything.

files = spm_select('FPListRec', directory, regularExpression);
assert(~isempty(files), 'Could not find any files in directory %s \n\n using regularExpression %s', directory, regularExpression)
files = cellstr(files);

%% Filter
% Filter our matches using the Filt string. Search the full paths of the
% matched files for the string "Filt". Only return matches that contain
% "Filt" somewhere within the full path.

boolFilt = contains(files, Filter);
if ~any(boolFilt) 
    disp(files)
    error('Filt %s filtered out all matches in files. Consider revising.', Filter)
end
files    = files(boolFilt);

end 