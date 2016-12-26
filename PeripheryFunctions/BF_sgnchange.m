function out = BF_sgnchange(y,doFind)
% BF_sgnchange      Where a data vector changes sign
%
%---INPUTS:
% y, the input vector
% doFind, (i) 0: returns a logical vector with 1s where the input changes sign
%         (ii) 1: returns a vector of indicies of where the input vector changes
%                 sign

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
    doFind = 0; % by default don't find, just return logical of where the sign changes (faster)
end

out = (y(2:end).*y(1:end-1) < 0); % successive differences change sign

% [case of adding equality, as <= 0, but I think better to
% make hard changes in sign as < 0]

if doFind
    out = find(out); % return indicies of where sign changes
end

end
