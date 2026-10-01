function out = CO_FirstCrossing(y, corrFun, threshold, whatOut)
% CO_FirstCrossing   The first crossing of the autocorrelation function across a threshold.
%
% Returns the lag at which the autocorrelation function (ACF) first crosses a
% threshold, as the first whole lag at which it is below the threshold (if it starts
% above) and as a point linearly interpolated between that lag and the lag before.
% Lags are in samples, counted from lag 0. If the ACF never crosses, both are N-1.
%
% ---INPUTS:
% y, the input time series
% corrFun, the self-correlation function to measure:
%          'ac': the linear autocorrelation function, calculated with CO_AutoCorr
%                (the default; the only option)
% threshold, the value to cross (default 0). Examples: 0 [first zero-crossing],
%            1/exp(1) [first 1/e crossing].
% whatOut, what to return (default 'both'): 'both' for a structure with the two
%          fields below, 'discrete' for the first, or 'continuous' for the second.
%
% ---OUTPUTS:
% firstCrossing, the first whole lag at which the ACF is on the other side of the
%       threshold,
% pointOfCrossing, the linearly interpolated lag of the crossing.

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
% Check inputs:
% ------------------------------------------------------------------------------

N = length(y); % the length of the time series

if nargin < 2 || isempty(corrFun)
	corrFun = 'ac'; % linear autocorrelation function
end
if nargin < 3 || isempty(threshold)
	threshold = 0;
end
if nargin < 4 || isempty(whatOut)
	whatOut = 'both';
end

% -------------------------------------------------------------------------------
% Select the self-correlation function as an inline function
% Eventually could add additional self-correlation functions
switch corrFun
	case 'ac'
		% Autocorrelation at all time lags
		% (a little inefficient to not do it incrementally, but maybe Fourier method)
		% (does it anyway...)
		corrs = CO_AutoCorr(y, [], 'Fourier');
	otherwise
		error('Unknown correlation function ''%s''', corrFun);
end

% -------------------------------------------------------------------------------
% Calculate point of crossing:
[firstCrossingIndex, pointOfCrossingIndex] = BF_PointOfCrossing(corrs, threshold);

% -------------------------------------------------------------------------------
% Assemble the appropriate output (structure or double):
% Convert from index space (1,2,…) to lag space (0,1,2,…):
switch whatOut
	case 'both'
		out.firstCrossing = firstCrossingIndex - 1;
		out.pointOfCrossing = pointOfCrossingIndex - 1;
	case 'discrete'
		out = firstCrossingIndex - 1;
	case 'continuous'
		out = pointOfCrossingIndex - 1;
end

end
