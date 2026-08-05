function out = DN_Moments(y, theMom, doNormalize)
% DN_Moments    A moment of the distribution of the input time series.
%
% ---INPUTS:
% y, the input data vector
% theMom, the moment to calculate (a scalar)
% doNormalize, whether to normalize by std(y)^theMom, giving the proper
%              scale-invariant standardized moment (true, default -- e.g.,
%              theMom=3 is skewness, theMom=4 is kurtosis), or to return
%              the raw, unnormalized central moment (false).
%
% Uses the moment function from Matlab's Statistics Toolbox.
%
% NOTE: prior to 2026-08, this always divided by std(y)^1 regardless of
% theMom, which is neither the raw central moment nor a scale-invariant
% standardized moment -- it has no statistical meaning beyond the special
% case where y is already unit-variance (where it happens to coincide
% with the standardized moment, since std(y)^1 = std(y)^theMom = 1). Any
% code relying on that specific (unintended) behavior on non-unit-variance
% input should now pass doNormalize=false and account for the change.

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
