function out = CO_AddNoise(y,tau,amiMethod,extraParam,randomSeed)
% CO_AddNoise
%
% Analyzes changes in the automutual information function with the addition of
% noise to the input time series.
%
% Adds Gaussian-distributed noise to the time series with increasing standard
% deviation, eta, across the range eta = 0, 0.1, ..., 2, and measures the
% mutual information at each point
% Can be measured using histograms with extraParam bins (implemented using
% CO_HistogramAMI), or using the Information Dynamics Toolkit.
%
% The output is a set of statistics on the resulting set of automutual
% information estimates, including a fit to an exponential decay, since the
% automutual information decreases with the added white noise.
%
% Can calculate these statistics for time delays 'tau', and for a number 'extraParam'
% bins.
%
% This algorithm is quite different, but was based on the idea of 'noise
% titration' presented in: "Titration of chaos with added noise", Chi-Sang Poon
% and Mauricio Barahona P. Natl. Acad. Sci. USA, 98(13) 7107 (2001)
%
%---INPUTS:
%
% y, the input time series
%
% tau, the time delay for computing AMI (using CO_HistogramAMI)
%
% amiMethod, the method for computing AMI (using CO_HistogramAMI)
%
% extraParam, e.g., the number of bins input to CO_HistogramAMI
%
% randomSeed: settings for resetting the random seed for reproducible results
%               (using BF_ResetSeed)

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
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
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Preliminary checks
% ------------------------------------------------------------------------------
% Check a curve-fitting toolbox license is available:
BF_CheckToolbox('curve_fitting_toolbox');

doPlot = 0; % plot outputs to figure

% ------------------------------------------------------------------------------
%% Check inputs
% ------------------------------------------------------------------------------
if nargin < 2
    tau = []; % set default in CO_HistogramAMI
end
if nargin < 3
    amiMethod = 'even'; % using evenly spaced bins in CO_HistogramAMI
end
if nargin < 4
    extraParam = []; % number of bins for CO_HistogramAMI
end
if nargin < 5
    randomSeed = [];
end

% ------------------------------------------------------------------------------
% Preliminaries
% ------------------------------------------------------------------------------
noiseRange = (0:0.1:2); % compare properties across this noise range
BF_ResetSeed(randomSeed); % reset the random seed if specified
numRepeats = length(noiseRange);
amis = zeros(numRepeats,1);
noise = randn(size(y));

% ------------------------------------------------------------------------------
% Compute the automutual information across a range of noise levels
% ------------------------------------------------------------------------------
% The *same* noise vector, noise, is added to the signal, with increasing
% standard deviation (one could imagine repeating the calculation with different
% random seeds)...
switch amiMethod
case {'std1','std2','quantiles','even'}
    % histogram-based methods using my naive implementation in CO_Histogram.m
    for i = 1:numRepeats
        amis(i) = CO_HistogramAMI(y+noiseRange(i)*noise,tau,amiMethod,extraParam);
    end
case {'gaussian','kernel','kraskov1','kraskov2'}
    for i = 1:numRepeats
        amis(i) = IN_AutoMutualInfo(y+noiseRange(i)*noise,tau,amiMethod,extraParam);
    end
end

% ------------------------------------------------------------------------------
% Statistics
% ------------------------------------------------------------------------------

% Proportion decreases:
out.pdec = sum(diff(amis) < 0)/(numRepeats-1);

% Mean change in AMI:
out.meanch = mean(diff(amis));

% Autocorrelation of AMIs:
out.ac1 = CO_AutoCorr(amis,1,'Fourier');
out.ac2 = CO_AutoCorr(amis,2,'Fourier');

% Noise level required to reduce ami to 1/e original value:
out.ecorrLength = noiseRange(find(amis<amis(1)/exp(1),1,'first'));

% First time amis drop under value x: 1, 0.5, 0.2, 0.1
firstUnderVals = [1,0.5,0.2,0.1];
for i = 1:length(firstUnderVals)
    out.(sprintf('firstUnder%u',firstUnderVals(i)*10)) = firstUnder_fn(firstUnderVals(i),noiseRange,amis);
end

% AMI at actual noise levels: 0.5, 1, 1.5 and 2
noiseLevels = [0.5,1,1.5,2];
for i = 1:length(noiseLevels)
    out.(sprintf('ami_at_%u',noiseLevels(i)*10)) = amis(noiseRange==noiseLevels(i));
end

% ------------------------------------------------------------------------------
% Fit exponential decay to output using Curve Fitting Toolbox
% ------------------------------------------------------------------------------
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[amis(1) -1]);
f = fittype('a*exp(b*x)','options',s);
[c, gof] = fit(noiseRange',amis,f);

% Output statistics on fit to an exponential decay
out.fitexpa = c.a;
out.fitexpb = c.b;
out.fitexpr2 = gof.rsquare;
out.fitexpadjr2 = gof.adjrsquare;
out.fitexprmse = gof.rmse;

% ------------------------------------------------------------------------------
% Fit linear function to output
% ------------------------------------------------------------------------------
p = polyfit(noiseRange',amis,1);
out.fitlina = p(1);
out.fitlinb = p(2);
linfit = polyval(p,noiseRange);
out.mse = mean((linfit' - amis).^2);

% ------------------------------------------------------------------------------
% Number of times the AMI function crosses its mean
% ------------------------------------------------------------------------------
out.pcrossmean = sum(BF_sgnchange(amis-mean(amis)))/(numRepeats-1);

% ------------------------------------------------------------------------------
% Plot output:
% ------------------------------------------------------------------------------
if doPlot
    figure('color','w'); box('on');
    cc = BF_getcmap('set1',2,1);
    % figure('color','w');
    hold on; box('on')
    plot(noiseRange,c.a*exp(c.b*noiseRange),'color',cc{2},'linewidth',2)
    plot(noiseRange,amis,'.-','color',cc{1})
    xlabel('\eta'); ylabel('AMI_1')
end


% ------------------------------------------------------------------------------
function firsti = firstUnder_fn(x,m,p)
    %% Find m for the first time p goes under x
    firsti = m(find(p < x,1,'first'));
    % If it never goes under -- saturate as NaN
    % (could alternatively set to the maximum m...)
    if isempty(firsti)
        firsti = NaN;
    end
end

end
