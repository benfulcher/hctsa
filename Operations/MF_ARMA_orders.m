function out = MF_ARMA_orders(y,pr,qr)
% MF_ARMA_orders    Compares a range of ARMA models fitted to a time series.
%
% Given a set of AR orders, p, and a set of MA orders, q, this operation fits
% ARMA(p,q) models to the time series and evaluates the goodness of fit from all
% combinations (p,q).
%
% Uses functions iddata, armax, and aic from Matlab's System Identification toolbox
%
%---INPUTS:
% y, the input time series
% pr, a vector specifying the range of AR model orders to analyze
% qr, a vector specifying the range of MA model orders to analyze
%
%---OUTPUTS: statistics on the appropriateness of different types of models,
% including the goodness of fit from the best model, and the optimal orders of
% fitted ARMA(p,q) models.
%
% ** Future improvements **
% (1) May want to quantify where a particular order starts to stand out, i.e.,
% may be quite sensitive to wanting p > 2, but may be quite indiscriminate
% when it comes to the MA order.
% (2) May want to do some prediction and get more statistics on quality of
% model rather than just in-sample FPE/AIC...

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
%% Check that a System Identification Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('identification_toolbox')

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
% Convert y to time series object
y = iddata(y,[],1);

% ARMA(p,q): p range, pr
if nargin < 2 || isempty(pr)
   pr = (1:10);
end

% ARMA(p,q): q range, qr
if nargin < 3 || isempty(qr)
    qr = (1:5);
end

% ------------------------------------------------------------------------------
%% Preliminaries
% ------------------------------------------------------------------------------
fpes = zeros(length(pr),length(qr));
aics = zeros(length(pr),length(qr));

% ------------------------------------------------------------------------------
%% Fit the models
% ------------------------------------------------------------------------------
for i = 1:length(pr)
    p = pr(i);
    for j = 1:length(qr)
        q = qr(j);

        % Fit the ARMA(p,q) model
        m = armax(y,[p,q]);

        % Get statistics on it
        fpes(i,j) = m.EstimationInfo.FPE;
        aics(i,j) = aic(m);
    end
end

% global minimum, aic
out.aic_min = min(aics(:));
[pi_opt, qi_opt] = find(aics == min(aics(:)),1,'first');
out.p_aic_opt = pr(pi_opt);
out.q_aic_opt = qr(qi_opt);

out.std_all_aics = std(aics(:)); % no idea why.
out.mean_all_aics = mean(aics(:)); % no idea why.

out.meanstd_aicsp = mean(std(aics));
out.meanstd_aicsq = mean(std(aics,[],2));


end
