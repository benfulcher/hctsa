function nlogL = DN_nlogL_norm(y)
% DN_nlogL_norm     Negative log likelihood of data coming from a Gaussian distribution.
%
% Fits a Gaussian distribution to the data using the normfit function in
% Matlab's Statistics Toolbox and returns the negative log likelihood of the
% data coming from a Gaussian distribution using the normlike function.
%
%---INPUT:
% y, a vector of data

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

[muhat, sigmahat] = normfit(y);
nlogL = normlike([muhat, sigmahat],y)/length(y);

% ** Somehow this just scales with length, regardless of the distribution of y?

end
