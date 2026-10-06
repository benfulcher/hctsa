function [ac1, ac2, runsz] = BF_ResidualStats(res, sstot)
% BF_ResidualStats   Statistics of the residuals of a fit, for remaining structure.
%
% Summarizes what a fitted model leaves unexplained, from the autocorrelation of the
% residuals at lags 1 and 2 (the 'Fourier' method of CO_AutoCorr) and a runs test on
% them (BF_RunsZ). The residuals are taken in the order given (time order, or order
% of increasing value of the fitted variable).
%
% A fit that is exact (the residual sum of squares is below 1e-12 of the total sum
% of squares) leaves only numerical error in the residuals, whose autocorrelation and
% runs are arbitrary and differ between implementations: all three outputs are then
% NaN instead of being computed from that noise.
%
% ---INPUTS:
% res, the residuals (column vector)
% sstot, the total sum of squares of the fitted data about its mean, which sets the
%       scale against which the residuals are judged to be negligible
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

if sum(res.^2) <= 1e-12 * sstot
	% exact fit: nothing but numerical error is left
	ac1 = NaN; ac2 = NaN; runsz = NaN;
	return
end

ac1 = CO_AutoCorr(res, 1, 'Fourier'); % autocorrelation of residuals at lag 1
ac2 = CO_AutoCorr(res, 2, 'Fourier'); % autocorrelation of residuals at lag 2
runsz = BF_RunsZ(res); % runs test on residuals

end
