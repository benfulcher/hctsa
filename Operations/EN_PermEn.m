% EN_PermEn
% 
% Computes the permutation entropy of order, ord, of a time series.
% 
% "Permutation Entropy: A Natural Complexity Measure for Time Series"
% C. Bandt and B. Pompe, Phys. Rev. Lett. 88(17) 174102 (2002)
% 
% Code is adapted from logisticPE.m code obtained from
% http://people.ece.cornell.edu/land/PROJECTS/Complexity/
% http://people.ece.cornell.edu/land/PROJECTS/Complexity/logisticPE.m
% 
% INPUTS:
% y, a time series
% ord, the order of permutation entropy
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
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

function out = EN_PermEn(y,ord)
% Ben Fulcher, 2009

if nargin < 2 || isempty(ord)
    ord = 2; % order 2
end

% Ensure y is a column vector
if size(y,1) > size(y,2);
    y = y';
end

% Use the Bruce Land and Damian Elias code to calculate the permutation entropy:
out = LA_permen(y,ord);

end