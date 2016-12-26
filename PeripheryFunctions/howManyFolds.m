function numFolds = howManyFolds(groupLabels,numClasses)
% Set the number of folds for k-fold cross validation using a heuristic
% (for small datasets with fewer than 10 examples per class):
%
%---INPUTS:
% groupLabels, labels of groups in the dataset
% numClasses, the number of classes of time series

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

if nargin < 2
    numClasses = max(groupLabels);
end

numPerClass = arrayfun(@(x)sum(groupLabels==x),1:numClasses);

% Make sure there are enough points in the smallest class to do proper
% cross-validation:
minPointsPerClass = min(numPerClass);

% Now the heuristic to set the number of folds:
if minPointsPerClass < 5
    numFolds = 2;
elseif minPointsPerClass < 10
    numFolds = 5;
else
    numFolds = 10;
end

end
