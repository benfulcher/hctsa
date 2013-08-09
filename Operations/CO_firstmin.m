% CO_FirstMin
% 
% Returns the time at which the first minimum in a given correlation function
% occurs.
% 
% INPUTS:
% y, the input time series
% minwhat, the type of correlation to minimize: either 'ac' for autocorrelation,
%           or 'mi' for automutual information
% 
% Note that selecting 'ac' is unusual operation: standard operations are the
% first zero-crossing of the autocorrelation (as in CO_FirstZero), or the first
% minimum of the mutual information function ('mi').
%
% The 'mi' option uses Rudy Moddemeijer's RM_information.m code that may or may
% not be great...
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
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

function out = CO_FirstMin(y,minwhat)
% Ben Fulcher, 2008

if nargin < 2 || isempty(minwhat)
    minwhat = 'mi'; % mutual information
end

N = length(y); % time-series length

switch minwhat
case 'mi'
    % automutual information implemented as RM_information
    corrfn = @(x) RM_information(y(1:end-x), y(1+x:end));
case 'ac'
    % autocorrelation implemented as CO_AutoCorr
    corrfn = @(x) CO_AutoCorr(y,x);
otherwise
    error('Unknown correlation function ''%s''',minwhat);
end

% Go through lags until find a minimum
a = zeros(N-1,1); % preallocate space for the autocorrelation function
for i = 0:N-1
    a(i+1) = corrfn(i); % calculate the auto-correlation at this lag
            
    if (i > 1) && (a(i-1) - a(i) > 0) && (a(i) - a(i+1) < 0); % minimum
        out = i-1; % I found the first minimum!
        return
    end
end

% If there is no minimum in the automutual information
out = N;

end

