function BF_CheckToolbox(theToolbox)
% BF_CheckToolbox   Check for a Matlab Toolbox license
%
% Checks for the given Matlab Toolbox license and produces an appropriate error
% message if there's a problem.
%
%---INPUT:
% theToolbox, the string identifying the toolbox (in format of license command
%               in Matlab)

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

% ------------------------------------------------------------------------------
% First get a more interpretable string for user feedback: theName
% ------------------------------------------------------------------------------
switch theToolbox
case 'statistics_toolbox'
    theName = 'Matlab''s Statistics Toolbox';

case 'identification_toolbox'
    theName = 'Matlab''s System Identification Toolbox';

case 'econometrics_toolbox'
    theName = 'Matlab''s Econometrics Toolbox';

case 'curve_fitting_toolbox'
    theName = 'Matlab''s Curve Fitting Toolbox';

case 'wavelet_toolbox'
    theName = 'Matlab''s Wavelet Toolbox';

case 'signal_toolbox'
    theName = 'Matlab''s Signal Processing Toolbox';

case 'database_toolbox'
    theName = 'Matlab''s Database Toolbox';

case 'robust_toolbox'
    theName = 'Matlab''s Robust Control Toolbox';

case 'financial_toolbox'
    theName = 'Matlab''s Financial Toolbox';

case 'distrib_computing_toolbox'
    theName = 'Matlab''s Parallel Computing Toolbox';

otherwise
    error('Unknown toolbox ''%s''\n',theToolbox);
end

% ------------------------------------------------------------------------------
% Now do the checks:
% ------------------------------------------------------------------------------
% 1. Check the toolbox exists in the current Matlab environment:
a = license('test',theToolbox);
if a == 0
    error('This function requires %s',theName);
end

% 2. Check to see if there's an available license for this toolbox:
[lic_free,~] = license('checkout',theToolbox); % Attempt to check out a license
if lic_free == 0
    error('Could not obtain a license for %s',theName);
end

end
