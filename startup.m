function startup()
% STARTUP   Add all paths required for the hctsa package.

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

% Resolve paths relative to this file's own location (not the current
% working directory), so startup works regardless of where it's called from:
hctsaDir = fileparts(mfilename('fullpath'));

% We use this function a bit:
addfcn = @(x) addpath(fullfile(hctsaDir,x));

fprintf(1,'Adding paths for the highly comparative time-series analysis package...\n')

% ------------------------------------------------------------------------------
%% First add all the basic paths:
% ------------------------------------------------------------------------------
addfcn('FeatureSets'); % files for setting input files corresponding to various feature sets
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
addpath(fullfile(hctsaDir,'Toolboxes','Danny_Kaplan'));

% Gaussian Process Toolbox, gpml, by Carl Edward Rasmussen & Hannes Nickisch:
fprintf(1,', Gaussian Process Code\n')
addpath(fullfile(hctsaDir,'Toolboxes','gpml'));
GP_startup % add nested directories

% Zoubin Gharamani's hmm toolbox, ZG_hmm
fprintf(1,'HMM toolbox')
addpath(fullfile(hctsaDir,'Toolboxes','ZG_hmm'));

% ARFIT Toolbox
fprintf(1,', ARfit toolbox')
addpath(fullfile(hctsaDir,'Toolboxes','ARFIT'));

% Michael Small's utilities
fprintf(1,', Michael Small\n')
addpath(fullfile(hctsaDir,'Toolboxes','Michael_Small'));

% DVV Toolbox
fprintf(1,'DVV Toolbox')
addpath(fullfile(hctsaDir,'Toolboxes','DVV_Toolbox'));

% Physionet
fprintf(1,', Physionet');
addpath(fullfile(hctsaDir,'Toolboxes','Physionet'));

% Max Little's steps/bumps toolbox
fprintf(1,', Max Little''s steps_bumps toolkit')
addpath(fullfile(hctsaDir,'Toolboxes','Max_Little','steps_bumps_toolkit'));

% Max Little's fastdfa code
fprintf(1,', fastdfa')
addpath(fullfile(hctsaDir,'Toolboxes','Max_Little','fastdfa'));

% Max Little's rpde code
fprintf(1,', rpde')
addpath(fullfile(hctsaDir,'Toolboxes','Max_Little','rpde'));

% nsamdf
fprintf(1,', nsamdf,\n');
addpath(fullfile(hctsaDir,'Toolboxes','nsamdf'));

% catch22
fprintf(1,'catch22')
addpath(fullfile(hctsaDir,'Toolboxes','catch22','wrap_Matlab'));

% Java information dynamics toolkit written by Joseph Lizier
% (should be ok to re-add this every time startup is run)
fprintf(1,', Information dynamics toolkit, ')
javaaddpath(fullfile(hctsaDir,'Toolboxes','infodynamics-dist','infodynamics.jar'));

% ------------------------------------------------------------------------------
% Add path for TISEAN binaries (compiled locally into Toolboxes/Tisean_3.0.1/bin
% by install.m; falls back to ~/bin for older manual installations):
% ------------------------------------------------------------------------------
tiseanBinaryLocation = fullfile(hctsaDir,'Toolboxes','Tisean_3.0.1','bin');
if ~isfolder(tiseanBinaryLocation)
    [~,homeDir] = system('echo $HOME'); % get system home directory
    homeDir = regexprep(homeDir,'[\s]',''); % remove whitespace
    tiseanBinaryLocation = fullfile(homeDir,'bin');
end
if isempty(regexp(getenv('PATH'),tiseanBinaryLocation,'once'))
    % Prepend (not append): an older manual install may already sit in
    % ~/bin (e.g., via a generic $HOME/bin PATH entry in .bash_profile
    % unrelated to TISEAN), which must not shadow the locally-built version.
    sysPath = [tiseanBinaryLocation,':',getenv('PATH')];
    setenv('PATH', sysPath)
    fprintf(1,'System path to TISEAN binaries: %s\n',tiseanBinaryLocation);
end

% ------------------------------------------------------------------------------
%% Finished:
% ------------------------------------------------------------------------------
fprintf(1,'\n---Done.\n')

end
