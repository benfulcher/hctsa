function out = TSTL_delaytime(y,maxDelay,past,randomSeed)
% TSTL_delaytime    Optimal delay time using the method of Parlitz and Wichard.
%
% Uses the TSTOOL code delaytime (this method is specified in the TSTOOL
% documentation but without reference).
%
% TSTOOL: http://www.physik3.gwdg.de/tstool/
%
%---INPUTS:
% y, column vector of time series data
%
% maxDelay, maximum value of the delay to consider (can also specify a
%           proportion of time series length)
%
% past, the TSTOOL documentation describes this parameter as "?", which is
%       relatively uninformative.
%
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%% Preliminaries
% ------------------------------------------------------------------------------
N = length(y); % length of time series
try
    s = signal(y); % convert to a signal for TSTOOL
catch
    error('Error converting time series to signal class using TSTOOL function ''signal''')
end

% ------------------------------------------------------------------------------
% Check Inputs:
% ------------------------------------------------------------------------------
% (1) Maximum delay, maxDelay
if nargin < 2 || isempty(maxDelay)
    maxDelay = 0.2; % 1/5 the length of the time series
end
if maxDelay < 1 && maxDelay > 0
    maxDelay = round(N*maxDelay); % specify a proportion of time series length
end

if maxDelay < 10,
    maxDelay = 10;
    fprintf(1,'Max delay set to its minimum: delaytime = 10\n');
end

% randomSeed: how to treat the randomization
if nargin < 4
    randomSeed = [];
end

% ------------------------------------------------------------------------------
%% Run
% ------------------------------------------------------------------------------

% Control the random seed (for reproducibility):
BF_ResetSeed(randomSeed);

% Run the TSTOOL delaytime function on the signal object time series:
tau = data(delaytime(s,maxDelay,past));

% plot(tau);

% ------------------------------------------------------------------------------
%% Output Statistics
% ------------------------------------------------------------------------------
% tau tends to start low and then rise to some (noisy) value
out.tau1 = tau(1);
out.tau2 = tau(2);
out.tau3 = tau(3);
out.difftau12 = tau(2)-tau(1);
out.difftau13 = tau(3)-tau(1);
out.meantau = mean(tau);
out.stdtau = std(tau);
out.mintau = min(tau);
out.maxtau = max(tau);

end
