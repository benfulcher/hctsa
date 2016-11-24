function out = DN_Fit_mle(y,fitWhat)
% DN_Fit_mle    Maximum likelihood distribution fit to data.
%
% Fits either a Gaussian, Uniform, or Geometric distribution to the data using
% maximum likelihood estimation via the Matlab function mle
% from the Statistics Toolbox.
%
%---INPUTS:
% y, the input data vector.
% fitWhat, the type of fit to do: 'gaussian', 'uniform', or 'geometric'.
%
%---OUTPUTS: parameters from the fit.

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

%-------------------------------------------------------------------------------
% Check Inputs:
%-------------------------------------------------------------------------------
if nargin < 2
    fitWhat = 'gaussian'; % fit a Gaussian by default
end

%-------------------------------------------------------------------------------
% Do the fitting:
%-------------------------------------------------------------------------------
switch fitWhat
case 'gaussian'
	phat = mle(y);
	out.mean = phat(1); % mean of Gaussian fit
    out.std = phat(2); % std of Gaussian fit

case 'uniform' % turns out to be shit
    phat = mle(y,'distribution','uniform');
    out.a = phat(1);
    out.b = phat(2);

case 'geometric'
    out = mle(y,'distribution','geometric'); % just a single output

otherwise
    error('Invalid fit specifier, ''%s''',fitWhat)
end

end
