function weights = InverseProbWeight(Labels)
% Returns inverse probability weights (used for unbalanced classification problems)
%
%---INPUT:
% Labels: integer class label assigned to each observation
%
%---OUTPUT:
% weights: inverse probability weight assigned to each observation

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
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

classNames = unique(Labels);
numClasses = length(classNames);
numObs = length(Labels);

weights = zeros(size(Labels));

for x = 1:numClasses
	isClass = (Labels == classNames(x));
	weights(isClass) = numObs/sum(isClass);
end

end
