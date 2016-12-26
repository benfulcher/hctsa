function out = MF_ResidualAnalysis(e)
% MF_ResidualAnalysis   Analysis of residuals from a model fit.
%
% Given an input residual time series residuals, e, this function returns a
% structure with fields corresponding to statistical tests on the residuals.
% These are motivated by a general expectation of model residuals to be
% uncorrelated.
%
%---INPUT:
% e, should be raw residuals as prediction minus data (e = yp - y) as a column
%       vector.

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
%% Preliminaries
% ------------------------------------------------------------------------------

% Check that a System Identification Toolbox license is available (for spa):
BF_CheckToolbox('identification_toolbox')

if size(e,2) > size(e,1)
    e = e'; % make sure residuals are a column vector
end
if all(e>0)
    warning('Very weird that all model residuals are positive...')
elseif all(e<0)
    warning('Very weird that all model residuals are negative...')
end
N = length(e);

% ------------------------------------------------------------------------------
%% Basic statiatics on residuals, then zscore
% ------------------------------------------------------------------------------
out.meane = mean(e);
out.meanabs = mean(abs(e));
out.rmse = sqrt(mean(e.^2));
out.stde = std(e);
out.mms = abs(mean(e))+abs(std(e));
out.maxonmean = max(e)/abs(mean(e));

if std(e) == 0
    e = zeros(length(e),1);
else
    e = zscore(e);
end

% ------------------------------------------------------------------------------
%% Identify any low-frequency trends in residuals
% ------------------------------------------------------------------------------
% Look for any low-frequency trends -- extract summaries from power
% spectrum.
g = spa(e); % smoothed power spectrum
% p = etfe(e); % periodogram
gf = g.frequency;
gS = g.Spectrumdata(:);

% Normalize them
% this is like normalizing the residuals to unit variance
gS = gS / (sum(gS)*(gf(2)-gf(1)));

% Look at proportion of power in fifths
b = round(linspace(0,length(gf),6));
out.p1_5 = sum(gS(b(1)+1:b(2)))*(gf(2)-gf(1));
out.p2_5 = sum(gS(b(2)+1:b(3)))*(gf(2)-gf(1));
out.p3_5 = sum(gS(b(3)+1:b(4)))*(gf(2)-gf(1));
out.p4_5 = sum(gS(b(4)+1:b(5)))*(gf(2)-gf(1));
out.p5_5 = sum(gS(b(5)+1:b(6)))*(gf(2)-gf(1));

% ------------------------------------------------------------------------------
%% Analyze autocorrelation in residuals
% ------------------------------------------------------------------------------
% See if there are any linear correlations in residuals.
% Also see if any of these are abnormally large (i.e., may be remnant
% autocorrelation at some level, or may be a characteristic shape in this
% function...)
% Will output both raw values and values scaled by sqrt(length), as is
% normal (within a constant).
maxLag = 25;

acs = CO_AutoCorr(e,1:maxLag,'Fourier'); % autocorrelations
sqrtN = sqrt(N);

% Output first three acfs
out.ac1 = acs(1);
out.ac2 = acs(2);
out.ac3 = acs(3);
out.ac1n = abs(acs(1))*sqrtN; % units of 1/sqrtN from zero
out.ac2n = abs(acs(2))*sqrtN; % units of 1/sqrtN from zero
out.ac3n = abs(acs(3))*sqrtN; % units of 1/sqrtN from zero

% Median normalized distance from zero
out.acmnd0 = median(abs(acs))*sqrtN;
out.acsnd0 = std(abs(acs))*sqrtN;
out.propbth = sum(abs(acs) < 2.6/sqrtN)/maxLag;

% First time to get below the significance threshold
out.ftbth = find(abs(acs) < 2.6/sqrtN,1,'first');
if isempty(out.ftbth)
    out.ftbth = maxLag+1;
end

% Durbin-Watson test statistic (like AC1)
out.dwts = sum((e(2:end)-e(1:end-1)).^2) / sum(e.^2);

% ------------------------------------------------------------------------------
%% Linear structure in residuals
% ------------------------------------------------------------------------------
% Fit a linear model and see if it picks up any structure.
% There's also a suggestion in 'resid' documentation to fit an arx model to
% the output of resid -- looks for correlations between inputs and
% outputs, perhaps?

% Fit a zero-mean AR process to residuals using the ARFIT package:
emsg = '';
try
    [~, Aest, ~, SBC, FPE] = ARFIT_arfit(e, 1, 10, 'sbc', 'zero');
catch emsg
end

if ~isempty(emsg)
    % (strcmp(emsg.message,'Time series too short.') || strcmp(emsg.message,'Matrix must be positive definite.'))
    warning('Error fitting AR model to residuals using ARFIT package: %s\n',emsg.message)
    out.popt = NaN; % Optimum order
    out.minsbc = NaN; % Best sbc
    out.minfpe = NaN; % Best fpe
    out.sbc1 = NaN; % SBC(1)
else
    out.popt = length(Aest); % Optimum order
    out.minsbc = min(SBC); % Best sbc
    out.minfpe = min(FPE); % Best fpe
    out.sbc1 = SBC(1);
end


% ------------------------------------------------------------------------------
%% Distribution tests
% ------------------------------------------------------------------------------
[~, p, ksstat] = kstest(e);
out.normksstat = ksstat;
out.normp = p;

end
