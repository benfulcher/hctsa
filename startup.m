addfcn = @(x) addpath(fullfile(pwd,x));
% addfcn = @(x) addpath(genpath(fullfile(pwd,x)));

fprintf(1,'Adding paths for the highly comparative time-series analysis package...')
addfcn('Database'); % code for setting up and communicating with the mySQL database
addfcn('Calculation'); % code for calculating results
addfcn('PlottingAnalysis'); % code for analysing and plotting results
addfcn('Operations'); % core code files for performing operations
addfcn('PeripheryFunctions'); % periphery functions used in the code toolbox
addfcn('TimeSeries'); % time series data files for analysis
fprintf(1,' Core directories added.\n')

%% Now add all Toolboxes:
fprintf(1,'Adding external time-series toolboxes...')
% Kaplan's routines:
fprintf(1,' Danny Kaplan')
addfcn('Toolboxes/Danny_Kaplan');

% Code by Marwan, from CRP Toolbox version 5.17  (R28.16)
fprintf(1,', Marwan')
addfcn('Toolboxes/Marwan_crptool');

% Gaussian Process Toolbox, gpml, by Carl Edward Rasmussen & Hannes Nickisch:
fprintf(1,', Gaussian Process Code')
addfcn('Toolboxes/gpml');
GP_startup % add nested directories

% ARFIT Toolbox
fprintf(1,', arfit tool')
addfcn('Toolboxes/arfit_tool');

% Michael Small's utilities
fprintf(1,', Michael Small')
addfcn('Toolboxes/Michael_Small');

% Zoubin Gharamani's hmm toolbox, ZG_hmm
fprintf(1,', HMM\n')
addfcn('Toolboxes/ZG_hmm');

% Code from Matlab Central
fprintf(1,'Matlab Central code')
addfcn('Toolboxes/MatlabCentral');

% Rudy Moddemeijer's code
fprintf(1,', Rudy Moddemeijer')
addfcn('Toolboxes/Rudy_Moddemeijer');

% Bill Davidson's hurst exponent code
fprintf(1,', Bill Davidson');
addfcn('Toolboxes/Bill_Davidson');

% Max Little's steps/bumps toolbox
fprintf(1,', steps_bumps toolkit')
addfcn('Toolboxes/Max_Little/steps_bumps_toolkit');

% Max Little's fastdfa code
fprintf(1,', fastdfa')
addfcn('Toolboxes/Max_Little/fastdfa');


% TSTOOL
fprintf(1,', TSTOOL')
addfcn('Toolboxes/OpenTSTOOL');
settspath('Toolboxes/OpenTSTOOL'); % this routine adds the necessary paths for OpenTSTOOL

fprintf(1,'. Done.\n')