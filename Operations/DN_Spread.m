function out = DN_Spread(y,spreadMeasure)
% DN_Spread     Measure of spread of the input time series.
%
% Returns the spread of the raw data vector, as the standard deviation,
% inter-quartile range, mean absolute deviation, or median absolute deviation.
%
%---INPUTS:
% y, the input data vector
%
% spreadMeasure, the spead measure:
%               (i) 'std': standard deviation
%               (ii) 'iqr': interquartile range
%               (iii) 'mad': mean absolute deviation
%               (iv) 'mead': median absolute deviation

% ------------------------------------------------------------------------------
% Copyright (C) 2017, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems (2017).
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
% Check Inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(spreadMeasure)
    spreadMeasure = 'std'; % return std by default
end

% ------------------------------------------------------------------------------
% Evaluate the spread measure
% ------------------------------------------------------------------------------
switch spreadMeasure
	case 'std'
        % Standard deviation
		out = std(y);

	case 'iqr'
        % Interquartile range
		out = iqr(y);

	case 'mad'
        % Mean absolute deviation
		out = mad(y,0);

    case 'mead'
        % Median absolute deviation
        out = mad(y,1);

    otherwise
        error('Unknown spread measure ''%s''',spreadMeasure)
end

end
