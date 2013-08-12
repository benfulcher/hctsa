% BF_CheckToolbox
% 
% Checks for the given Matlab Toolbox license and produces an appropriate error
% message if there's a problem.
% 
% INPUT:
% thetoolbox, the string identifying the toolbox (in format of license command
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

function BF_CheckToolbox(thetoolbox)

% First get a more interpretable string for user feedback: thename
switch thetoolbox
case 'identification_toolbox'
    thename = 'Matlab''s System Identification Toolbox';
    
case 'econometrics_toolbox'
    thename = 'Matlab''s Econometrics Toolbox';
    
case 'curve_fitting_toolbox'
    thename = 'Matlab''s Curve Fitting Toolbox';
    
case 'wavelet_toolbox'
    thename = 'Matlab''s Wavelet Toolbox';
    
case 'signal_toolbox'
    thename = 'Matlab''s Signal Processing Toolbox';
    
otherwise
    fprintf(1,'Unknown toolbox ''%s''\n',thetoolbox);
end

% Now do the checks:
% 1. Check the toolbox exists in the current Matlab environment:
a = license('test',thetoolbox);
if a == 0
    error('This function requires %s',thename);
end
    
% 2. Check to see if there's an available license for this toolbox:
[lic_free,~] = license('checkout',thetoolbox); % Attempt to check out a license
if lic_free == 0
    error('Could not obtain a license for %s');
end

end 