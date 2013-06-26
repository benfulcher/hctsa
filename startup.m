addfcn = @(x) addpath(fullfile(pwd,x));
% addfcn = @(x) addpath(genpath(fullfile(pwd,x)));

addfcn('Database');
addfcn('Operations');
addfcn('PeripheryFunctions');
addfcn('TimeSeries');
addfcn('TSQCoreRoutines');

%% Now add all Toolboxes:

% Kaplan's routines:
addfcn('Toolboxes/Kaplan');

% CRP Toolbox by Marwan, version 5.13, Release 26
addfcn('Toolboxes/crptool');

% TSTOOL
addfcn('Toolboxes/OpenTSTOOL');
settspath('Toolboxes/OpenTSTOOL'); % this routine adds the necessary paths for OpenTSTOOL

% Gaussian Processes
addfcn('Toolboxes/gpml');

% arfit Toolbox
addfcn('Toolboxes/arfit_tool');

% Michael Small's utilities
addfcn('Toolboxes/MSmall_utilities');

% Zoubin Gharamani's hmm toolbox
addfcn('Toolboxes/hmm');

% Max Little's steps/bumps toolbox
addfcn('Toolboxes/steps_bumps_toolkit');