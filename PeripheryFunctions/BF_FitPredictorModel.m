function [A, K, C] = BF_FitPredictorModel(y, model, order)
% BF_FitPredictorModel   Fit a linear time-series model and return it in innovations form.
%
% Fits the model used by MF_steps_ahead and MF_CompareTestSets to the series y and returns
% it as x(t+1) = A x(t) + K e(t), y(t) = C x(t) + e(t), for predictions with BF_SSPredict.
%
%---INPUTS:
% y, the time series
% model, 'ar': an AR model fitted by forward-backward least squares (BF_ARfb), with
%               order an integer or 'best' (the order from 1 to 10 minimizing Schwarz's
%               Bayesian criterion, from ARFIT_arfit with a zero mean);
%        'arma': an ARMA model fitted by conditional least squares (BF_ARMAfit),
%               with order = [p, q];
%        'ss': a state-space model identified by N4SID (BF_SubspaceID), with order an
%               integer or 'best'.
% order, the model order (see above)
%
% Errors if the model cannot be fitted.

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

y = y(:);
switch model
    case 'ar'
        if ischar(order) && strcmp(order, 'best')
            [~, Aest] = ARFIT_arfit(y, 1, 10, 'sbc', 'zero');
            order = length(Aest);
        end
        a = BF_ARfb(y, order);
        [A, K, C] = BF_ARMAtoSS(a, 1);
    case 'arma'
        [a, c] = BF_ARMAfit(y, order(1), order(2));
        [A, K, C] = BF_ARMAtoSS(a, c);
    case 'ss'
        m = BF_SubspaceID(y, order);
        A = m.A; K = m.K; C = m.C;
    otherwise
        error('Unknown model ''%s''', model);
end

end
