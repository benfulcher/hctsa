function out = SY_StdNthDer(y, n)
% SY_StdNthDer   Standard deviation of the nth derivative of the time series.
%
% Based on an idea by Vladimir Vassilevsky, a DSP and Mixed Signal Design
% Consultant, in a comp.soft-sys.matlab (MATLAB newsgroup) posting, who stated that
% "You can measure the standard deviation of the nth derivative, if you like".
%
% The derivative is estimated very simply by taking successive increments of the
% time series (diff(y, n), with a unit time step); the process is repeated to
% obtain higher order derivatives.
%
% Note that this idea is popular in the heart-rate variability literature, cf.
% M. Brennan, M. Palaniswami and P. Kamen, "Do existing measures of Poincare plot
% geometry reflect nonlinear features of heart rate variability?", IEEE Trans. Biomed.
% Eng. 48(11), 1342-1347 (2001) (and function MD_hrv_classic). DOI: 10.1109/10.959330
%
% ---INPUTS:
% y, the time series to analyze
%
% n, the order of derivative to analyze (default: 2)
%
% ---OUTPUTS:
% a scalar: the standard deviation of the nth difference of the time series

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

if nargin < 2 || isempty(n)
	n = 2;
end

yd = diff(y, n); % crude method of taking a derivative that could be improved
% upon in future

if isempty(yd)
	error('Time series (N = %u) too short to compute differences at %u', ...
		  length(y), n);
end
out = std(yd);

end
