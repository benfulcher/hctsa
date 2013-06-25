% addfcn = @(x) addpath(genpath(fullfile(pwd,x)));
addpath('Database');
addpath('Operations');
addpath('PeripheryFunctions');
addpath('TimeSeries');
addpath('TSQCoreRoutines');

%% Now add all Toolboxes:

% Kaplan's routines:
addpath('Toolboxes/Kaplan');

% CRP Toolbox by Marwan, version 5.13, Release 26
addpath('Toolboxes/crptool');

% TSTOOL
addpath('Toolboxes/OpenTSTOOL');
settspath('Toolboxes/OpenTSTOOL'); % this routine adds the necessary paths for OpenTSTOOL

% Gaussian Processes
addpath('Toolboxes/gpml');

% arfit Toolbox
addpath('Toolboxes/arfit_tool');

% Michael Small's utilities
addpath('Toolboxes/MSmall_utilities');

% Zoubin Gharamani's hmm toolbox
addpath('Toolboxes/hmm');

% Max Little's steps/bumps toolbox
addpath('Toolboxes/steps_bumps_toolkit');