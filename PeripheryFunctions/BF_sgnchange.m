% BF_sgnchange
% 
% Returns where the input vector, y, changes sign.
% 
% INPUTS:
% y, the input vector
% dofind, (i) 0: returns a logical vector with 1s where the input changes sign
%         (ii) 1: returns a vector of indicies of where the input vector changes
%                 sign
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

function out = BF_sgnchange(y,dofind)

if nargin < 2
    dofind = 0; % by default don't find, just return logical of where the sign changes (faster)
end

out = (y(2:end).*y(1:end-1) < 0); % successive differences change sign

% [case of adding equality, as <= 0, but I think better to 
% make hard changes in sign as < 0]

if dofind
    out = find(out); % return indicies of where sign changes
end

end