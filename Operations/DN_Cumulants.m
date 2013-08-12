% DN_Cumulants
% 
% Very simple function that uses the skewness and kurtosis functions in 
% Matlab's Statistics Toolbox to calculate these higher order moments of input time series, y
% 
% INPUTS:
% 
% y, the input time series
% 
% whatcum, the type of higher order moment:
%           (i) 'skew1', skewness
%           (ii) 'skew2', skewness correcting for bias
%           (iii) 'kurt1', kurtosis
%           (iv) 'kurt2', kurtosis correcting for bias
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

function out = DN_Cumulants(y,whatcum)
% Ben Fulcher, 2008

if nargin < 2 || isempty(whatcum)
    whatcum = 'skew1'; % do skewness by default
end

switch whatcum
case 'skew1' % skewness
	out = skewness(y);
    
case 'skew2' % corrects for bias
    out = skewness(y,0);
    
case 'kurt1' % kurtosis
	out = kurtosis(y);
    
case 'kurt2' % corrects for bias
    out = kurtosis(y,0);        
    
otherwise
    error('Unknown cumulant ''%s'' specified: should be ''skew1'', ''skew2'', ''kurt1'', or ''kur2''',whatcum)
end

end