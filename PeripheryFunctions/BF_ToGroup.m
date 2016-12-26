function group = BF_ToGroup(groupIndices,maxLength)
% BF_ToGroup    convert from cell form to vector form of group labels
%
% Converts to a vector of group labels from a cell, where each element is a
% vector of indices for that group (or vice-versa).

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

if iscell(groupIndices) % Convert to vector of group indices

    if all(cellfun(@isempty,groupIndices)) % all empty -- rubbish
        error('Nothing in the group indices')
    end

    % Second input only important for cell -> vector transformation
    if nargin < 2 || isempty(maxLength)
        maxLength = max(cellfun(@max,groupIndices));
        % Make the length of output equal to the maximum index
    end

    group = zeros(maxLength,1);

    for i = 1:length(groupIndices)
        group(groupIndices{i}) = i;
    end

    if sum(cellfun(@length,groupIndices)) < length(group)
        warning('Group is missing some labels')
    end

else
    % Convert from vector of group labels to cell of group indices:

    numGroups = max(groupIndices);
    group = cell(numGroups,1);
    for i = 1:numGroups
        group{i} = find(groupIndices==i);
    end
end

end
