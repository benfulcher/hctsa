% ------------------------------------------------------------------------------
% startup.m
% ------------------------------------------------------------------------------
% Add all paths required for the highly comparative time-series analysis package

% We use this function a bit:
addfcn = @(x) addpath(fullfile(pwd,x));

fprintf(1,'Adding paths for the highly comparative time-series analysis package...')

% ------------------------------------------------------------------------------
%% First add all the basic paths:
% ------------------------------------------------------------------------------
addfcn('Database'); % code for setting up and communicating with the mySQL database
addfcn('Calculation'); % code for calculating results
addfcn('PlottingAnalysis'); % code for analysing and plotting results
addfcn('Operations'); % core code files for performing operations
addfcn('PeripheryFunctions'); % periphery functions used in the code toolbox
addfcn('TimeSeries'); % time series data files for analysis
fprintf(1,' Core directories added.\n')

% ------------------------------------------------------------------------------
%% Now add all the external code packages and periphery toolboxes:
% ------------------------------------------------------------------------------
fprintf(1,'Adding external time-series toolboxes...')
% Kaplan's routines:
fprintf(1,' Danny Kaplan')
addpath(fullfile(pwd,'Toolboxes','Danny_Kaplan'));

% Code by Marwan, from CRP Toolbox version 5.17  (R28.16)
fprintf(1,', Marwan')
addpath(fullfile(pwd,'Toolboxes','Marwan_crptool'));

% Gaussian Process Toolbox, gpml, by Carl Edward Rasmussen & Hannes Nickisch:
fprintf(1,', Gaussian Process Code')
addpath(fullfile(pwd,'Toolboxes','gpml'));
GP_startup % add nested directories

% Zoubin Gharamani's hmm toolbox, ZG_hmm
fprintf(1,', HMM toolbox\n')
addpath(fullfile(pwd,'Toolboxes','ZG_hmm'));

% ARFIT Toolbox
fprintf(1,'ARfit toolbox')
addpath(fullfile(pwd,'Toolboxes','ARFIT'));

% Michael Small's utilities
fprintf(1,', Michael Small')
addpath(fullfile(pwd,'Toolboxes','Michael_Small'));

% Code from Matlab Central
fprintf(1,', Matlab Central code')
addpath(fullfile(pwd,'Toolboxes','MatlabCentral'));

% Rudy Moddemeijer's code
fprintf(1,', Rudy Moddemeijer')
addpath(fullfile(pwd,'Toolboxes','Rudy_Moddemeijer'));

% Land and Elias (code from http://people.ece.cornell.edu/land/PROJECTS/Complexity/)
fprintf(1,', Land and Elias');
addpath(fullfile(pwd,'Toolboxes','Land_and_Elias'));

% TS Research
fprintf(1,', TS_Research\n')
addpath(fullfile(pwd,'Toolboxes','TS_Research'));

% Physionet
fprintf(1,', Physionet');
addpath(fullfile(pwd,'Toolboxes','Physionet'));

% Max Little's steps/bumps toolbox
fprintf(1,', steps_bumps toolkit')
addpath(fullfile(pwd,'Toolboxes','Max_Little','steps_bumps_toolkit'));

% Max Little's fastdfa code
fprintf(1,', fastdfa')
addpath(fullfile(pwd,'Toolboxes','Max_Little','fastdfa'));

% TSTOOL
fprintf(1,', TSTOOL')
addpath(fullfile(pwd,'Toolboxes','OpenTSTOOL'));
% Run the routine adds the necessary paths for OpenTSTOOL:
settspath(fullfile(pwd,'Toolboxes','OpenTSTOOL'));

% ------------------------------------------------------------------------------
% Add path for TISEAN (ASSUMING in ~/bin directory):
% ------------------------------------------------------------------------------
[~,homeDir] = system('echo $HOME'); % get system home directory
homeDir = regexprep(homeDir,'[\s]',''); % remove white space
tiseanBinaryLocation = fullfile(homeDir,'bin');
if isempty(regexp(getenv('PATH'),tiseanBinaryLocation))
    sysPath = [getenv('PATH'),':',tiseanBinaryLocation];
    setenv('PATH', sysPath)
    fprintf(1,',\nSystem path to TISEAN binaries: %s',tiseanBinaryLocation);
    clear sysPath
end
clear homeDir tiseanBinaryLocation

% TISEAN also requires this DYLD path to be set (I assume this works also on Windows):
setenv('DYLD_LIBRARY_PATH','/usr/local/bin');
fprintf(1,', DYLD library path');

% ------------------------------------------------------------------------------
%% Finished:
% ------------------------------------------------------------------------------
clear addfcn % clear the add function
fprintf(1,'\n---Done.\n')