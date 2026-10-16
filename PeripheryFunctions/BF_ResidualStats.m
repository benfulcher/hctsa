function [ac1, ac2, runsz] = BF_ResidualStats(res, y)
% BF_ResidualStats   Statistics of the residuals of a fit, for remaining structure.
%
% Summarizes what a fitted model leaves unexplained, from the autocorrelation of the
% residuals at lags 1 and 2 (the 'Fourier' method of CO_AutoCorr) and a runs test on
% them (BF_RunsZ). The residuals are taken in the order given (time order, or order
% of increasing value of the fitted variable).
%
% All three outputs are NaN if the residuals have no variation beyond floating-point
% round-off, as for an exact fit: the autocorrelation and runs of round-off are not
% meaningful (and differ between implementations). The criterion is a standard
% deviation of the residuals of at most 10*N*eps times the root-mean-square of the
% fitted data (N values), a bound on the rounding error of a fit to them: a fitted
% value accumulates rounding over sums of up to N terms (at most about N*eps/2 of the
% magnitude of the data), and a sinusoid evaluated at times up to N has a phase
% rounding error of up to about 1.5*pi*N*eps (see PP_Compare). Small but genuine
% residuals, such as the quantization noise left by an exact fit to 24-bit data, are
% computed.
%
% ---INPUTS:
% res, the residuals (column vector)
% y, the fitted data (column vector), whose magnitude sets the round-off scale
%
% ---OUTPUTS:
% ac1, ac2, the autocorrelation of the residuals at lags 1 and 2
% runsz, the signed z-statistic of a runs test on the residuals (BF_RunsZ)

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

if std(res) <= 10 * length(res) * eps * sqrt(mean(y.^2))
	% nothing but round-off is left (an exact fit)
	ac1 = NaN; ac2 = NaN; runsz = NaN;
	return
end

ac1 = CO_AutoCorr(res, 1, 'Fourier'); % autocorrelation of residuals at lag 1
ac2 = CO_AutoCorr(res, 2, 'Fourier'); % autocorrelation of residuals at lag 2
runsz = BF_RunsZ(res); % runs test on residuals

end
