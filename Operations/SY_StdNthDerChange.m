function out = SY_StdNthDerChange(y,maxd)
% SY_StdNthDerChange    How the output of SY_StdNthDer changes with order parameter.
%
% Order parameter controls the derivative of the signal.
%
% Operation inspired by a comment on the Matlab Central forum: "You can
% measure the standard deviation of the n-th derivative, if you like." --
% Vladimir Vassilevsky, DSP and Mixed Signal Design Consultant from
% http://www.mathworks.de/matlabcentral/newsreader/view_thread/136539
%
%---INPUTS:
% y, the input time series
%
% maxd, the maximum derivative to take.
%
%---OUTPUTS:
% An exponential function, f(x) = Aexp(bx), is fitted to the variation across
% successive derivatives; outputs are the parameters and quality of this fit.
%
% Typically an excellent fit to exponential: regular signals decrease, irregular
% signals increase...?

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

% ------------------------------------------------------------------------------
%% Check that a Curve-Fitting Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('curve_fitting_toolbox')

doPlot = 0; % plot outputs

if nargin < 2 || isempty(maxd)
    maxd = 10; % do 10 by default
end

ms = zeros(maxd,1);
for i = 1:maxd
    ms(i) = SY_StdNthDer(y,i);
end

if doPlot
    figure('color','w'); box('on');
    plot(ms,'o-k')
end

% Fit exponential growth/decay using the Curve-Fitting Toolbox
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1, 0.5*sign(ms(end)-ms(1))]);
f = fittype('a*exp(b*x)','options',s);
[c, gof] = fit((1:maxd)',ms,f);
out.fexp_a = c.a;
out.fexp_b = c.b; % this is important
out.fexp_r2 = gof.rsquare; % this is more important!
out.fexp_adjr2 = gof.adjrsquare;
out.fexp_rmse = gof.rmse;

end
