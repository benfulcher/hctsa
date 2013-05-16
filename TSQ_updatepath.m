% Update Path
% Sets all the path information for all necessary routines, functions, metrics, data,...
disp('<<>>TSQ_updatepath<<>>: adding necessary paths for data, routines, and functions...');

% startbit = '/Volumes/fulcher/MATLAB/Time_Series_Resources/'; % MAC
startbit = 'H:\MATLAB\Time_Series_Resources\';

% Directories in which data files are kept
% (can help to have the data split between multiple directories so that
% 	no individual directory has more than a few thousand files, which slow file accesses)
addpath([startbit 'Data_1']);
addpath([startbit 'Data_2']);
addpath([startbit 'Data_3']);

% Add core routines
addpath([startbit 'TSQ_core_routines\']);

% Metric function files
addpath([startbit 'OperationFiles']);

% Subroutines for use in general code
addpath([startbit 'PeripheryFunctions']);

% % Matlab functions and routines for data analysis
% addpath([startbit 'TSQ_routines\ClusteringClassification']);
% addpath([startbit 'TSQ_routines\DataMining']);
% addpath([startbit 'TSQ_routines\Housekeeping']);
% addpath([startbit 'TSQ_routines\Medical']);
% addpath([startbit 'TSQ_routines\Visualization']);
% addpath([startbit 'TSQ_routines\Workshop']);
% addpath([startbit 'TSQ_routines\Tests']);
% addpath([startbit 'midlists']);

%% Other toolboxes
% CRP Toolbox by Marwan, version 5.13, Release 26
addpath([startbit 'Toolboxes\crptool']);

% TSTOOL
addpath([startbit 'Toolboxes\OpenTSTOOL']);
settspath([startbit 'Toolboxes\OpenTSTOOL']); % this routine adds the necessary paths for OpenTSTOOL

% Gaussian Processes
addpath([startbit 'Toolboxes\gpml']);

% arfit Toolbox
addpath([startbit 'Toolboxes\arfit_tool']);

% Michael Small's utilities
addpath([startbit 'Toolboxes\MSmall_utilities']);

% Zoubin Gharamani's hmm, discrete hmm toolboxes
addpath([startbit 'Toolboxes\hmm']);
% addpath([startbit 'Toolboxes\dhmm']);

% Max Little's steps/bumps toolbox
addpath([startbit 'Toolboxes\steps_bumps_toolkit']);

% Voice Analysis Toolkit from Thanasis
% addpath([startbit 'Toolboxes\Voice_analysis_toolkit']);

% Matlab implementation of LASSO, LARS, the elastic net and SPCA
% http://www2.imm.dtu.dk/pubdb/views/publication_details.php?id=3897
addpath([startbit 'Toolboxes\LASSO']);

% Isomap nonlinear dimensionality reduction
addpath([startbit 'Toolboxes\IsoMap']);

% LibSVM
% addpath([startbit 'Toolboxes\LibSVM']);

% DSPCA
% addpath([startbit 'Toolboxes\DSPCA']);

% FastICA
% addpath([startbit 'Toolboxes\FastICA']);

%% Spider
% if addspider
% startbit = 'H:\MATLAB\Time_Series_Resources\';
% addpath([startbit 'spider'])
% use_spider;
% end
