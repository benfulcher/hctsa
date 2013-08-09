% CO_TSTL_AutoCorrMethod
% 
% Estimates the autocorrelation function using a fast Fourier Transform method
% implemented in TSTOOL and returns the mean square discrepancy between the
% autocorrelation coefficients obtained in this way from those obtained in the
% time domain using CO_AutoCorr.
% 
% TSTOOL: http://www.physik3.gwdg.de/tstool/
% 
% No real rationale behind this, other than the difference in autocorrelations
% computed by the two methods may somehow be informative of something about the
% time series...
% 
% INPUTS:
% y, the input time series
% maxlag, the maximum time lag to compute up to -- will compare autocorrelations
%         up to this value
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

function out = CO_TSTL_AutoCorrMethod(y,maxlag)
% Ben Fulcher, October 2009

if nargin < 2 || isempty(maxlag)
    maxlag = 50; % compare across the first maxlag autocorrelations
end

% First maxlag autocorrelations
co_fft = data(acf(signal(y),maxlag*2));

nlags = length(co_fft);
co_ben = zeros(nlags,1);
for i = 1:nlags
    co_ben(i) = CO_AutoCorr(y,i-1);
end

out = norm(co_ben - co_fft);

end