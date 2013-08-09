% DN_CustomSkewness
% 
% Estimates custom skewness measures, the Pearson and Bowley skewnesses.
% 
% INPUTS:
% y, the input time series
% 
% whichskew, the skewness measure to calculate, either 'pearson' or 'bowley'
% 
% The Bowley skewness uses the quantile function from Matlab's Statistics Toolbox
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

function out = DN_CustomSkewness(y,whichskew)
% Ben Fulcher, 2009

switch whichskew
    case 'pearson'
        out = (3*mean(y)-median(y))./std(y);
        
    case 'bowley'
        qs = quantile(y,[0.25, 0.5, 0.75]);
        out = (qs(3)+qs(1) - 2*qs(2))./(qs(3)-qs(1));
        % Quartile skewness coefficient
        
    otherwise
        error('Unknown skewness type ''%s''',whichskew)
end

end