function out = EN_DistributionEntropy(y,histOrKS,numBins,olremp)
% EN_DistributionEntropy    Distributional entropy.
%
% Estimates of entropy from the distribution of a data vector. The
% distribution is estimated either using a histogram with numBins bins, or as a
% kernel-smoothed distribution, using the ksdensity function from Matlab's
% Statistics Toolbox with width parameter, w (specified as the iunput numBins).
%
% An optional additional parameter can be used to remove a proportion of the
% most extreme positive and negative deviations from the mean as an initial
% pre-processing.
%
%---INPUTS:
%
% y, the input time series
%
% histOrKS: 'hist' for histogram, or 'ks' for ksdensity
%
% numBins: (*) (for 'hist'): an integer, uses a histogram with that many bins
%          (*) (for 'ks'): a positive real number, for the width parameter for
%                       ksdensity (can also be empty for default width
%                                       parameter, optimum for Gaussian)
%
% olremp [opt]: the proportion of outliers at both extremes to remove
%               (e.g., if olremp = 0.01; keeps only the middle 98% of data; 0
%               keeps all data. This parameter ought to be less than 0.5, which
%               keeps none of the data).
%               If olremp is specified, returns the difference in entropy from
%               removing the outliers.

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

doPlot = 0; % plot outputs to figure

% ------------------------------------------------------------------------------
%% Check inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(histOrKS)
    histOrKS = 'hist'; % use histogram by default
end
if nargin < 3 % (can be empty for default width for ksdensity)
    numBins = 10; % use 10 bins
end
if nargin < 4
    olremp = 0;
end

% ------------------------------------------------------------------------------
% (1) Remove outliers?
% ------------------------------------------------------------------------------
if olremp ~= 0
    yHat = y(y >= quantile(y,olremp) & y <= quantile(y,1-olremp));
    if isempty(yHat)
        % removed the entire time series?!
        % shouldn't be possible for good values of olremp with equality
        % in the above inequalities
        out = NaN; return
    else
        % Return the difference in entropy from removing outliers
        out = EN_DistributionEntropy(y,histOrKS,numBins) - ...
                EN_DistributionEntropy(yHat,histOrKS,numBins);
        return
    end
end

% ------------------------------------------------------------------------------
% (2) Form the histogram
% ------------------------------------------------------------------------------
switch histOrKS
case 'hist' % Use histogram to calculate pdf
    if isnumeric(numBins)
        [px,binEdges] = histcounts(y,numBins,'Normalization','probability');
    else
        [px,binEdges] = histcounts(y,'BinMethod',numBins,'Normalization','probability');
    end
    % Compute bin centers:
    xr = mean([binEdges(1:end-1); binEdges(2:end)]);
    % Compute bin widths:
    binWidths = diff(binEdges);

case 'ks' % Use ksdensity to calculate pdf
    if isempty(numBins)
        [px, xr] = ksdensity(y,'function','pdf'); % selects optimal width
    else
        [px, xr] = ksdensity(y,'width',numBins,'function','pdf'); % uses specified width
    end
    binWidths = ones(1,length(px))*(xr(2)-xr(1));

otherwise
    error('Unknown distribution method -- specify ''ks'' or ''hist''') % error; must specify 'ks' or 'hist'
end

if doPlot
    figure('color','w'); box('on');
    plot(xr,px,'k')
end

% ------------------------------------------------------------------------------
% (3) Compute the entropy sum and return it as output
% ------------------------------------------------------------------------------
% 0*log0 = 0:
out = -sum(px(px>0).*log(px(px>0)./binWidths(px>0)));

end
