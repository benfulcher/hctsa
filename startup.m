addfcn = @(x) addpath(fullfile(pwd,x));
% addfcn = @(x) addpath(genpath(fullfile(pwd,x)));

fprintf(1,'Adding paths for the highly comparative time-series analysis package...')
addfcn('Database');
addfcn('Operations');
addfcn('PeripheryFunctions');
addfcn('TimeSeries');
addfcn('TSQCoreRoutines');
fprintf(1,' core directories added...\n')

%% Now add all Toolboxes:
fprintf(1,'Adding external time-series toolboxes...')

% Kaplan's routines:
fprintf(1,' Kaplan')
addfcn('Toolboxes/Kaplan');

% CRP Toolbox by Marwan, version 5.13, Release 26
fprintf(1,', crptool')
addfcn('Toolboxes/crptool');

% Gaussian Processes
fprintf(1,', gpml')
addfcn('Toolboxes/gpml');

% arfit Toolbox
fprintf(1,', arfit')
addfcn('Toolboxes/arfit_tool');

% Michael Small's utilities
fprintf(1,', Michael Small')
addfcn('Toolboxes/MSmall_utilities');

% Zoubin Gharamani's hmm toolbox
fprintf(1,', HMM')
addfcn('Toolboxes/hmm');

% Max Little's steps/bumps toolbox
fprintf(1,', steps_bumps')
addfcn('Toolboxes/steps_bumps_toolkit');

% TSTOOL
fprintf(1,', TSTOOL')
addfcn('Toolboxes/OpenTSTOOL');
settspath('Toolboxes/OpenTSTOOL'); % this routine adds the necessary paths for OpenTSTOOL

fprintf(1,'. Done.\n')