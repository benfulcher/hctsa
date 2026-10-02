function out = DN_PLeft(y, th)
% DN_PLeft   Distance from the mean beyond which a given proportion of the data lie.
%
% Finds the distance from the mean that is exceeded by a proportion th of the
% data points, and normalizes it by the standard deviation of the data. That
% is, quantile(abs(y - mean(y)), 1 - th)/std(y), using the quantile function
% from MATLAB's Statistics Toolbox. Deviations above and below the mean are
% pooled (it could be generalized to treat them separately). For a Gaussian
% distribution and th = 0.05, the output is about 1.96.
%
% ---INPUTS:
% y, the input data vector
% th, the proportion of data points further from the mean than the output
%       distance (default: 0.1)
%
% ---OUTPUTS:
% a scalar: the distance from the mean exceeded by a proportion th of the
% data, in units of the standard deviation.

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

if nargin < 2 || isempty(th)
	th = 0.1; % default
end

p = quantile(abs(y - mean(y)), 1 - th);

% A proportion, th, of the data lie further than p from the mean
out = p / std(y);

end
