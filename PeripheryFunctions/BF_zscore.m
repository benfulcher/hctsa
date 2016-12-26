function zscoredData = BF_zscore(inputData)
% BF_zscore     z-score the input data vector (without a Statistics Toolbox licence)
%
%---INPUT:
% inputData, the input time series (or any vector).
%
%---OUTPUT:
% zscoredData, the z-scored transformation of the input.

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

% Check for NaNs:
if any(isnan(inputData))
    error('inputData contains NaNs');
end

% By default, z-score twice to reduce the numerical error:
zscoredData = (inputData - mean(inputData)) / std(inputData);
zscoredData = (zscoredData - mean(zscoredData)) / std(zscoredData);

end
