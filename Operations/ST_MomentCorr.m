function out = ST_MomentCorr(x, windowLength, wOverlap, mom1, mom2, whatTransform)
% ST_MomentCorr   Correlations between simple statistics in local windows of a time series.
%
% Slides a window along the series and computes two statistics (mom1 and mom2) in
% each window, then measures how related they are across windows: their
% correlation, and the density of windows in the plane of the two statistics. The
% series can first be transformed (whatTransform).
%
% The idea to implement this was that of Prof. Nick S. Jones (Imperial College
% London).
%
% ---INPUTS:
% x, the input time series
%
% windowLength, the sliding window length, in samples (a value below 1 is taken as
%       a proportion of the time-series length, rounded up; default: 0.02)
%
% wOverlap, the overlap between consecutive windows, in samples (a value below 1
%       is taken as a fraction of the window length, rounded down; default: 1/5)
%
% mom1, mom2: the statistics to investigate correlations between (in each window):
%               (i) 'iqr': interquartile range
%               (ii) 'median': median
%               (iii) 'std': standard deviation (about the local mean)
%               (iv) 'mean': mean
%           (defaults: mom1 = 'mean', mom2 = 'std')
%
% whatTransform, the transformation to apply to the time series before analyzing it:
%               (i) 'abs': takes absolute values of all data points
%               (ii) 'sqrt': takes the square root of absolute values of all
%                            data points
%               (iii) 'sq': takes the square of every data point
%               (iv) 'none': does no transformation (default)
%
% ---OUTPUTS:
% R, the correlation coefficient between mom1 and mom2 across windows
% absR, the absolute value of R
% density, the density of windows in the (mom1, mom2) plane: the number of windows
%       divided by the area of the box bounding them, range(mom1)*range(mom2)

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

doPlot = false; % plot outputs

N = length(x); % number of samples in the input signal

% -------------------------------------------------------------------------------
% Check inputs, set defaults
% -------------------------------------------------------------------------------

% Sliding window length (samples)
if nargin < 2 || isempty(windowLength)
	windowLength = 0.02; % 2% of the time-series length
end
if windowLength < 1
	windowLength = ceil(N * windowLength);
end

% Sliding window overlap length
if nargin < 3 || isempty(wOverlap)
	wOverlap = 1 / 5;
end
if wOverlap < 1 % specify a fraction OF THE WINDOW LENGTH
	wOverlap = floor(windowLength * wOverlap);
end

if nargin < 4 || isempty(mom1)
	mom1 = 'mean';
end
if nargin < 5 || isempty(mom2)
	mom2 = 'std';
end
if nargin < 6 || isempty(whatTransform)
	whatTransform = 'none';
end

% ------------------------------------------------------------------------------
% Apply the specified whatTransformormation:
% ------------------------------------------------------------------------------
switch whatTransform
	case 'abs'
		x = abs(x);
	case 'sq'
		x = x.^2;
	case 'sqrt'
		x = sqrt(abs(x));
	case 'none'
		% x = x;
	otherwise
		error('Unknown tranformation ''%s''', whatTransform)
end

% ------------------------------------------------------------------------------
% Create the windows:
% ------------------------------------------------------------------------------
x_buff = buffer(x, windowLength, wOverlap, 'nodelay');
numWindows = (N / (windowLength - wOverlap)); % number of windows

if size(x_buff, 2) > numWindows
	% fprintf(1,'Should have %u columns but we have %u: removing last one',numWindows,size(x_buff,2))
	x_buff = x_buff(:, 1:end - 1); % lose last point
end
pointsPerWindow = size(x_buff, 1);
if pointsPerWindow == 1
	error('This time series (N = %u) is too short to extract %u windows.', N, numWindows);
end

% ok, now we have the sliding window ('buffered') signal, x_buff
% first calculate the first moment in all the windows (each column is a
% 'window' of the signal
M1 = SUB_CalcMeMoments(x_buff, mom1);
M2 = SUB_CalcMeMoments(x_buff, mom2);

R = corrcoef(M1, M2);
out.R = R(2, 1); % correlation coefficient
out.absR = abs(R(2, 1)); % absolute value of correlation coefficient
out.density = length(M1) / (range(M1) * range(M2)); % density of points in M1--M2 space: (number of windows) / (bounding-box area)

if doPlot
	figure('color', 'w');
	plot(M1, M2, '.k');
end

% ------------------------------------------------------------------------------
function moms = SUB_CalcMeMoments(x_buff, momType)
	switch momType
		case 'mean'
			moms = mean(x_buff);
		case 'std'
			moms = std(x_buff);
		case 'median'
			moms = median(x_buff);
		case 'iqr'
			moms = iqr(x_buff);
		otherwise
			error('Unknown statistic ''%s''.', momType)
	end
end
% ------------------------------------------------------------------------------

end
