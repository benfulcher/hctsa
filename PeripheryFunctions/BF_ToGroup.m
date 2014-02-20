% ------------------------------------------------------------------------------
% BF_ToGroup
% ------------------------------------------------------------------------------
%
% Converts to Matlab's group form (a vector of group labels) from my
% GroupIndices form (which has indicies for each group as a cell).
% 
% Or vice-versa (6/1/2014)
% 
% --------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function group = BF_ToGroup(GroupIndices,MaxLength)

if iscell(GroupIndices) % Convert to vector of group indices
    
    if all(cellfun(@isempty,GroupIndices)) % all empty -- rubbish
        error('Nothing in the group indices')
    end
    
    % Second input only important for cell -> vector transformation
    if nargin < 2 || isempty(MaxLength)
        MaxLength = max(cellfun(@max,GroupIndices));
        % Make the length of output equal to the maximum index
    end
    
    group = zeros(MaxLength,1);

    for i = 1:length(GroupIndices)
        group(GroupIndices{i}) = i;
    end

    if sum(cellfun(@length,GroupIndices)) < length(group)
        warning('Group is missing some labels')
    end
    
else
    % Convert from vector of group labels to cell of group indices:
    
    NGroups = max(GroupIndices);
    group = cell(NGroups,1);
    for i = 1:NGroups
        group{i} = find(GroupIndices==i);
    end
end

end