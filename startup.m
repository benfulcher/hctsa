% STARTUP   Add all paths required for the hctsa package.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

% We use this function a bit:
addfcn = @(x) addpath(fullfile(pwd,x));

fprintf(1,'Adding paths for the highly comparative time-series analysis package...\n')

% ------------------------------------------------------------------------------
%% First add all the basic paths:
% ------------------------------------------------------------------------------
addfcn('Database'); % code for setting up and communicating with the mySQL database
addfcn('Calculation'); % code for calculating results
addfcn('PlottingAnalysis'); % code for analysing and plotting results
addfcn('Operations'); % core code files for performing operations
addfcn('PeripheryFunctions'); % periphery functions used in the code toolbox
addfcn('TimeSeries'); % time series data files for analysis
fprintf(1,'Core directories added.\n')

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
fprintf(1,', Gaussian Process Code\n')
addpath(fullfile(pwd,'Toolboxes','gpml'));
GP_startup % add nested directories

% Zoubin Gharamani's hmm toolbox, ZG_hmm
fprintf(1,'HMM toolbox')
addpath(fullfile(pwd,'Toolboxes','ZG_hmm'));

% ARFIT Toolbox
fprintf(1,', ARfit toolbox')
addpath(fullfile(pwd,'Toolboxes','ARFIT'));

% Michael Small's utilities
fprintf(1,', Michael Small')
addpath(fullfile(pwd,'Toolboxes','Michael_Small'));

% Code from Matlab Central
fprintf(1,', Matlab Central code')
addpath(fullfile(pwd,'Toolboxes','MatlabCentral'));

% Rudy Moddemeijer's code
fprintf(1,', Rudy Moddemeijer\n')
addpath(fullfile(pwd,'Toolboxes','Rudy_Moddemeijer'));

% DVV Toolbox
fprintf(1,'DVV Toolbox')
addpath(fullfile(pwd,'Toolboxes','DVV_Toolbox'));

% Physionet
fprintf(1,', Physionet');
addpath(fullfile(pwd,'Toolboxes','Physionet'));

% Max Little's steps/bumps toolbox
fprintf(1,', Max Little''s steps_bumps toolkit')
addpath(fullfile(pwd,'Toolboxes','Max_Little','steps_bumps_toolkit'));

% Max Little's fastdfa code
fprintf(1,', fastdfa')
addpath(fullfile(pwd,'Toolboxes','Max_Little','fastdfa'));

% Max Little's rpde code
fprintf(1,', rpde')
addpath(fullfile(pwd,'Toolboxes','Max_Little','rpde'));

% Misc code
fprintf(1,', misc')
addpath(fullfile(pwd,'Toolboxes','Misc'));

% TSTOOL
fprintf(1,', TSTOOL\n')
addpath(fullfile(pwd,'Toolboxes','OpenTSTOOL'));
% Run the routine adds the necessary paths for OpenTSTOOL:
settspath(fullfile(pwd,'Toolboxes','OpenTSTOOL'));

% Java information dynamics toolkit written by Joseph Lizier
% (should be ok to re-add this every time startup is run)
fprintf(1,'Information dynamics toolkit, ')
javaaddpath(fullfile(pwd,'Toolboxes','infodynamics-dist','infodynamics.jar'));

% ------------------------------------------------------------------------------
% Add path for TISEAN (ASSUMING in ~/bin directory):
% ------------------------------------------------------------------------------
[~,homeDir] = system('echo $HOME'); % get system home directory
homeDir = regexprep(homeDir,'[\s]',''); % remove whitespace
tiseanBinaryLocation = fullfile(homeDir,'bin');
if isempty(regexp(getenv('PATH'),tiseanBinaryLocation))
    sysPath = [getenv('PATH'),':',tiseanBinaryLocation];
    setenv('PATH', sysPath)
    fprintf(1,'System path to TISEAN binaries: %s\n',tiseanBinaryLocation);
    clear sysPath
end
clear homeDir tiseanBinaryLocation

% TISEAN also requires this DYLD path to be set (I assume this works also on
% cygwin/Windows?):
setenv('DYLD_LIBRARY_PATH','/usr/local/bin');
fprintf(1,'DYLD library path (for TISEAN).');

% ------------------------------------------------------------------------------
%% Finished:
% ------------------------------------------------------------------------------
clear addfcn % clear the add function
fprintf(1,'\n---Done.\n')
