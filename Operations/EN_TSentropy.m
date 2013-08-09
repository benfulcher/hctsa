% TSentropy
% 
% Estimates the Tsallis entropy of a signal using a parameter q that
% measures the non-extensivity of the system; q = 1 recovers the Shannon
% entropy.
% 
% INPUTS:
% x, the time series
% q, the non-extensivity parameter
% 
% Uses code written by D. Tolstonogov and obtained from
% http://download.tsresearchgroup.com/all/tsmatlablink/TSentropy.m.
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

function out = EN_TSentropy(x, q)
% Wrapper for TS_entropy
% Ben Fulcher 2009

if nargin < 2
    q = 1; % Shannon entropy by default
end

out = TS_entropy(x,q);

end