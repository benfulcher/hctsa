function out = SY_StdNthDerChange(y, maxd)
% SY_StdNthDerChange   How the output of SY_StdNthDer changes with the order of the derivative.
%
% Computes SY_StdNthDer(y, n) (the standard deviation of the nth difference of the
% time series) for orders n = 1, ..., maxd, and characterizes how it varies with n
% in two ways: by an exponential fit, and directly through the order at which it
% is smallest.
%
% Operation inspired by a comment on the MATLAB Central forum: "You can measure the
% standard deviation of the n-th derivative, if you like." (Vladimir Vassilevsky,
% DSP and Mixed Signal Design Consultant), from
% http://www.mathworks.de/matlabcentral/newsreader/view_thread/136539
%
% An exponential function, f(x) = a*exp(b*x), is fitted to the variation across
% successive derivatives: regular signals decrease, irregular signals increase.
% This exponential-decay/growth picture only holds when std(diff(y,n)) is monotonic
% across n. Many real (especially oversampled/smooth) series instead show
% successive differencing REDUCE std up to some order (removing trend or
% nonstationary drift) before over-differencing increases it again: a classic
% Box-Jenkins ARIMA-order-selection U-shape that a monotonic exponential cannot
% represent (on a 20-series sample of the Bonn EEG dataset, 20/20 showed this
% interior minimum, with a median exponential-fit r^2 of only 0.11). The minOrder,
% minOrderInterp, minRatio, overDiffRatio, and isInterior outputs characterize this
% directly, alongside the exponential fit. Needs the Curve Fitting Toolbox; if the
% exponential fit fails, the fexp_* outputs are NaN and the others are still returned.
%
% ---INPUTS:
% y, the input time series
%
% maxd, the maximum derivative (difference) order to take (default: 10)
%
% ---OUTPUTS:
% fexp_a, fexp_b, fexp_r2, fexp_rmse: the parameters a and b, the R^2, and the
%       root-mean-square error of the exponential fit f(n) = a*exp(b*n)
% minOrder, the order (1 to maxd) at which the standard deviation is smallest
% minOrderInterp, that order refined between integers by a parabola through the
%       three points around the minimum (equal to minOrder if the minimum is at
%       either end)
% minRatio, the smallest standard deviation divided by that at order 1
% overDiffRatio, the standard deviation at order maxd divided by the smallest
% isInterior, 1 if the minimum is strictly between order 1 and maxd (a U shape),
%       otherwise 0

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
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

% Check that a Curve-Fitting Toolbox license is available:
BF_CheckToolbox('curve_fitting_toolbox');

doPlot = false; % plot outputs

% -------------------------------------------------------------------------------
% Set defaults:
% -------------------------------------------------------------------------------
if nargin < 2 || isempty(maxd)
	maxd = 10; % do 10 by default
end
% -------------------------------------------------------------------------------

ms = zeros(maxd, 1);
for i = 1:maxd
	ms(i) = SY_StdNthDer(y, i);
end

if doPlot
	figure('color', 'w'); box('on');
	plot(ms, 'o-k')
end

% Fit exponential growth/decay using the Curve-Fitting Toolbox
s = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [1, 0.5 * sign(ms(end) - ms(1))]);
f = fittype('a*exp(b*x)', 'options', s);
try
	[c, gof] = fit((1:maxd)', ms, f);
	out.fexp_a = c.a;
	out.fexp_b = c.b; % this is important
	out.fexp_r2 = gof.rsquare; % this is more important!
	% fexp_adjr2 dropped: near-duplicate of fexp_r2 (adjustment is for 2 free
	% parameters against maxd=10 data points, so it barely moves).
	out.fexp_rmse = gof.rmse;
catch
	% The fit failed (e.g., non-finite values): NaN for the fit fields, but
	% still report the directly computed minimum-order statistics below
	out.fexp_a = NaN;
	out.fexp_b = NaN;
	out.fexp_r2 = NaN;
	out.fexp_rmse = NaN;
end

% ------------------------------------------------------------------------------
%% Directly characterize the minimum-variance differencing order
% ------------------------------------------------------------------------------
% (complements the exponential fit above, which can't represent a U-shaped
% curve -- see NOTE in the header)
[minStd, minInd] = min(ms);
out.minOrder = minInd;
out.minRatio = minStd / ms(1); % how much differencing helped, relative to order 1
out.overDiffRatio = ms(end) / minStd; % how much std rises again past the optimum
out.isInterior = double(minInd > 1 && minInd < maxd); % genuine U-shape vs. monotonic

% Sub-integer refinement of the minimizing order via parabolic interpolation
% of the three points straddling the discrete minimum -- a cheap, smarter
% alternative to repeating the sweep with fractional-order differencing:
if out.isInterior
	y0 = ms(minInd - 1); y1 = ms(minInd); y2 = ms(minInd + 1);
	denom = y0 - 2 * y1 + y2;
	if denom ~= 0
		out.minOrderInterp = minInd + 0.5 * (y0 - y2) / denom;
	else
		out.minOrderInterp = minInd;
	end
else
	out.minOrderInterp = minInd;
end

end
