% ST_FitPolynomial
% 
% Fits a polynomial of order k to the time series, and returns the mean
% square error of the fit.
% 
% Usually kind of a stupid thing to do with a time series, but it's sometimes
% somehow informative for time series with large trends.
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

function out = ST_FitPolynomial(y,k)
% Ben Fulcher, 2009

N = length(y); % the length of the time series (number of samples)
t = (1:N)'; % Get a range for the time axis for time series y

% Supress the (valid!) warning from stupidly fitting a polynomial to a time series...
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
cf = polyfit(t,y,k);
warning('on','MATLAB:polyfit:RepeatedPointsOrRescale');

f = polyval(cf,t);
out = sum((y-f).^2)/N; % mean RMS ERROR OF FIT

% % n=10;
% errs=zeros(n,1);
% x=1:length(y);
% for i=1:n
% cf=polyfit(x,y',i);
% f=polyval(cf,x);
% errs(i)=sum((y'-f).^2);
% % end
% % hold on;plot(f,'k')
% % plot(errs);

end