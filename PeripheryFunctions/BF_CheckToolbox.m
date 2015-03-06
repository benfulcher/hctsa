% ------------------------------------------------------------------------------
% BF_CheckToolbox
% ------------------------------------------------------------------------------
% 
% Checks for the given Matlab Toolbox license and produces an appropriate error
% message if there's a problem.
% 
%---INPUT:
% theToolbox, the string identifying the toolbox (in format of license command
%               in Matlab)
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function BF_CheckToolbox(theToolbox)

% ------------------------------------------------------------------------------
% First get a more interpretable string for user feedback: theName
% ------------------------------------------------------------------------------
switch theToolbox
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