% BF_label_space
% 
% Returns a cell that only keeps every 'skip' objects from the input string, and
% makes all other entries empty.
% Useful for labeling plots with too many elements to be legible otherwise.
% 
% INPUTS:
% k, the input cell of strings
% skip, the increment at which string labels will be kept
% 
% OUTPUT:
% labels, the transformed labels
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function labels = BF_label_space(k,skip)
% Ben Fulcher, 2008, 2009

if nargin < 2
    skip = 1; % skip not specified: don't skip any
end 
l = length(k);

labels = cell(l,1);
for i = 1:l
    if mod(i,skip) == 0
        labels{i} = k{i};
    else
        labels{i} = '';
    end
end

end