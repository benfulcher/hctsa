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
%% (the product name that ver reports), and the toolbox's folder: theFolder
%-------------------------------------------------------------------------------
switch theToolbox
case 'statistics_toolbox'
    theName = 'Statistics and Machine Learning Toolbox';
    theFolder = 'stats';
case 'curve_fitting_toolbox'
    theName = 'Curve Fitting Toolbox';
    theFolder = 'curvefit';
case 'signal_toolbox'
    theName = 'Signal Processing Toolbox';
    theFolder = 'signal';
case 'identification_toolbox'
    theName = 'System Identification Toolbox';
    theFolder = 'ident';
case 'wavelet_toolbox'
    theName = 'Wavelet Toolbox';
    theFolder = 'wavelet';
case 'econometrics_toolbox'
    theName = 'Econometrics Toolbox';
    theFolder = 'econ';
case 'robust_toolbox'
    theName = 'Robust Control Toolbox';
    theFolder = 'robust';
case 'financial_toolbox'
    theName = 'Financial Toolbox';
    theFolder = 'finance';
case 'database_toolbox'
    theName = 'Database Toolbox';
    theFolder = 'database';
case 'parallel_computing_toolbox'
    theName = 'Parallel Computing Toolbox';
    theFolder = 'parallel';
case 'distrib_computing_toolbox'
    theName = 'Parallel Computing Toolbox';
    theFolder = 'parallel';
otherwise
    error('Unknown toolbox ''%s''.\n',theToolbox);
end

%-------------------------------------------------------------------------------
%% Now do the checks:
%-------------------------------------------------------------------------------
% 1. Check the toolbox exists in the current Matlab environment:
%
% Installation is checked with ver(theFolder), which reads the toolbox's
% Contents.m from the Matlab path, matching the product name. It is not checked
% with matlab.addons.installedAddons, which can return an empty or incomplete
% list in a newly started Matlab session: on a fresh 48-worker parallel pool
% (NCI Gadi, R2025a), 44/48 workers got an empty list 2s after starting while
% ver, exist and license all answered correctly, and an incomplete list (once
% cached) failed every later call for the missing toolboxes on that worker.
%
% Both checks are read-only queries whose answer can't change during a session
% (unlike license('checkout',...) below, which reserves a license seat -- a real
% side effect, so that one is left uncached), so each toolbox's result is cached
% here: this function is called by ~42 operations.
persistent cachedKeys cachedResults
if isempty(cachedKeys)
    cachedKeys = {};
    cachedResults = false(0);
end
if doInstallCheck
    cacheKey = ['installed:',theToolbox];
else
    cacheKey = ['licensed:',theToolbox];
end

outFlag = false;
idx = find(strcmp(cacheKey,cachedKeys),1);
if isempty(idx)
    if doInstallCheck
        % have toolbox installed
        v = ver(theFolder);
        haveToolbox = any(strcmp({v.Name},theName));
    else
        % have toolbox license
        haveToolbox = logical(license('test',theToolbox));
    end
    cachedKeys{end+1} = cacheKey;
    cachedResults(end+1) = haveToolbox;
else
    haveToolbox = cachedResults(idx);
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
        if doInstallCheck
            error('This function requires %s but you don''t have it installed.',theName);
        else
            error('This function requires %s but no license for it was found.',theName);
        end
    end

    % 2. Check to see if there's an available license for this toolbox:
    [licenseFree,~] = license('checkout',theToolbox); % Attempt to check out a license
    if ~licenseFree
        error('This function requires %s but I could not obtain a license for it.',theName);
    end
end

end
