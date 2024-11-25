function [outFlag,theName] = BF_CheckToolbox(theToolbox,infoMode,doInstallCheck)
% BF_CheckToolbox   Check for a Matlab Toolbox license
%
% Checks for the given Matlab Toolbox license and produces an appropriate error
% message if there's a problem.
%
%---INPUTS:
%
% theToolbox, the string identifying the toolbox (in format of license command
%               in Matlab)
% infoMode, set true to just output information about a toolbox (not errors if missing).
% doInstallCheck, set true to check on installation (not license availability).
%
%---OUTPUTS:
%
% outFlag, true if the toolbox is NOT installed
% theName, text of the toolbox name

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

%-------------------------------------------------------------------------------
%% Check inputs
%-------------------------------------------------------------------------------
% By default, flag error if no toolbox
if nargin < 2 || isempty(infoMode)
    infoMode = false;
end
if nargin < 3 || isempty(doInstallCheck)
    doInstallCheck = true;
end

%-------------------------------------------------------------------------------
%% Get a more interpretable string for user feedback: theName
%% (and for testing against text in matlab.addons.installedAddons)
%-------------------------------------------------------------------------------
switch theToolbox
case 'statistics_toolbox'
    theName = 'Statistics and Machine Learning Toolbox';
case 'curve_fitting_toolbox'
    theName = 'Curve Fitting Toolbox';
case 'signal_toolbox'
    theName = 'Signal Processing Toolbox';
case 'identification_toolbox'
    theName = 'System Identification Toolbox';
case 'wavelet_toolbox'
    theName = 'Wavelet Toolbox';
case 'econometrics_toolbox'
    theName = 'Econometrics Toolbox';
case 'robust_toolbox'
    theName = 'Robust Control Toolbox';
case 'financial_toolbox'
    theName = 'Financial Toolbox';
case 'database_toolbox'
    theName = 'Database Toolbox';
case 'parallel_computing_toolbox'
    theName = 'Parallel Computing Toolbox';
case 'distrib_computing_toolbox'
    theName = 'Parallel Computing Toolbox';
otherwise
    error('Unknown toolbox ''%s''.\n',theToolbox);
end

%-------------------------------------------------------------------------------
%% Now do the checks:
%-------------------------------------------------------------------------------
% 1. Check the toolbox exists in the current Matlab environment:
outFlag = false;
if doInstallCheck
    installedAddOns = matlab.addons.installedAddons;
    % have toolbox installed
    haveToolbox = any(ismember(installedAddOns.Name,theName));
else
    % have toolbox license
    haveToolbox = license('test',theToolbox);
end
if infoMode
    % Just checking availability for info (e.g., during installation)
    if ~haveToolbox
        warning(['Some hctsa features require Matlab''s %s',...
                        ' but no installation was found.'],theName)
        outFlag = true;
    end
else
    % Want to use (and checkout) a license
    if ~haveToolbox
        error('This function requires %s but you don''t have it installed'.',theName);
    end

    % 2. Check to see if there's an available license for this toolbox:
    [licenseFree,~] = license('checkout',theToolbox); % Attempt to check out a license
    if ~licenseFree
        error('This function requires %s but I could not obtain a license for it).',theName);
    end
end

end
