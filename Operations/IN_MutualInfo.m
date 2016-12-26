function out = IN_MutualInfo(y1,y2,estMethod,extraParam)
% IN_MutualInfo     Mutual information of two data vectors.
%
% Uses the information dynamics toolkit implementation.
%
%---INPUTS:
%
% y1: input time series 1
% y2: input time series 2
%
% estMethod: the estimation method used to compute the mutual information:
%           (*) 'gaussian'
%           (*) 'kernel'
%           (*) 'kraskov1'
%           (*) 'kraskov2'
%
% cf. Kraskov, A., Stoegbauer, H., Grassberger, P., Estimating mutual
% information: http://dx.doi.org/10.1103/PhysRevE.69.066138

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
% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(y1) || isempty(y2)
    error('Need to provide two input vectors');
end
if size(y1)~=size(y2)
    error('Input vectors need to be the same size... :/');
end

if nargin < 3 || isempty(estMethod)
    estMethod = 'kernel';
end

if nargin < 4
    extraParam = [];
end

% ------------------------------------------------------------------------------
% Initialize miCalc object (don't add noise!):
miCalc = IN_Initialize_MI(estMethod,extraParam,0);

% Set observations to two time series:
miCalc.setObservations(y1, y2);

% Compute mutual information:
out = miCalc.computeAverageLocalOfObservations();

end
