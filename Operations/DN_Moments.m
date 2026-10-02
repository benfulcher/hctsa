function out = DN_Moments(y, theMom, doNormalize)
% DN_Moments   A moment of the distribution of the values of a time series.
%
% The central moment of order theMom of the values, ignoring their order in
% time, using the moment function from MATLAB's Statistics Toolbox. By
% default it is standardized by std(y)^theMom to give a scale-invariant
% quantity: theMom = 3 is then the skewness and theMom = 4 the kurtosis (3 for
% a Gaussian).
%
% ---INPUTS:
% y, the input data vector
% theMom, the order of the moment to calculate (a scalar)
% doNormalize, whether to normalize by std(y)^theMom, giving the
%       scale-invariant standardized moment (true, the default), or to return
%       the raw, unnormalized central moment (false)
%
% ---OUTPUTS:
% a scalar: the standardized or raw central moment of order theMom.
%
% ---NOTES:
% Earlier versions always divided by std(y)^1 regardless of theMom,
% which is neither the raw central moment nor a scale-invariant standardized
% moment: it has no statistical meaning beyond the special case where y is
% already unit-variance (where it coincides with the standardized moment,
% since std(y)^1 = std(y)^theMom = 1). Any code relying on that specific
% (unintended) behavior on non-unit-variance input should now pass
% doNormalize = false and account for the change.
%
% The moment function averages over N values, but std normalizes by N - 1, so
% for short series the standardized moments are slightly smaller than the
% textbook ones.

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

if nargin < 3 || isempty(doNormalize)
    doNormalize = true;
end

if doNormalize
    out = moment(y, theMom) / std(y)^theMom; % scale-invariant standardized moment
else
    out = moment(y, theMom); % raw central moment
end

end
