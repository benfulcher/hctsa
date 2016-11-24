function out = CO_TSTL_amutual(y,maxTau,numBins,versionTwo)
% CO_TSTL_amutual   Automutual information calculation using code from TSTOOL.
%
% Uses amutual code from TSTOOL, which uses a histogram method with n bins to
% estimate the mutual information of a time series across a range of
% time-delays, tau.
%
% TSTOOL: http://www.physik3.gwdg.de/tstool/
%
%---INPUTS:
%
% y, the time series
%
% maxTau, the maximum lag for which to calculate the auto mutual information
%
% numBins, the number of bins for histogram calculation
%
% versionTwo, uses amutual2 instead of amutual (from the TSTOOL package)
%
%---OUTPUTS:
% A number of statistics of the function over the range of tau, including the
% mean mutual information, its standard deviation, first minimum, proportion of
% extrema, and measures of periodicity in the positions of local maxima.

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
N = length(y); % length of time series
s = signal(y); % convert to signal object for TSTOOL

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(maxTau)
    maxTau = ceil(N/4);
end
maxTau0 = maxTau; % keep this original number if changed
maxTau = min(maxTau,ceil(N/2)); % Don't go above N/2

if nargin < 3 || isempty(numBins)
    numBins = round(sqrt(N/10)); % this is an arbitrary choice (!!) ;-)
end

if nargin < 4 || isempty(versionTwo)
    versionTwo = 0; % use amutual by default
end

% ------------------------------------------------------------------------------
%% Run
% ------------------------------------------------------------------------------
if versionTwo
    % Check existence of code:
    if ~exist('amutual2','file')
        error('''amutual2'' not found -- ensure the TSTOOL package is installed correctly??\n');
    end
    ami = data(amutual2(s,maxTau));
else
    % Check existence of code:
    if ~exist('amutual','file')
        error('''amutual'' not found -- ensure the TSTOOL package is installed correctly??\n');
    end
    ami = data(amutual(s,maxTau,numBins));
end

% Plot results
if doPlot
    figure('color','w'); box('on');
    plot(ami,'-ok');
end

% ------------------------------------------------------------------------------
% Change automutual information vector to a structure for output
% ------------------------------------------------------------------------------
for i = 1:maxTau0+1
    if i <= maxTau
        out.(sprintf('ami%u',i)) = ami(i);
    else
        out.(sprintf('ami%u',i)) = NaN;
    end
end

% ------------------------------------------------------------------------------
% Compute outputs
% ------------------------------------------------------------------------------
lami = length(ami);

% Mean mutual information over this lag range
out.mami = mean(ami);
out.stdami = std(ami);

% First miniimum of mutual information across range
dami = diff(ami);
extremai = find(dami(1:end-1).*dami(2:end) < 0);
out.pextrema = length(extremai)/(lami-1);
if isempty(extremai)
   out.fmmi = lami; % actually represents lag, because indexes don't but diff delays by 1
else
    out.fmmi = extremai(1);
end

% Periodicities in local maxima?
maximai = find(dami(1:end-1) > 0 & dami(2:end) < 0) + 1;
dmaximai = diff(maximai);
% is there a big peak in dmaxima?
% (no need to normalize since a given method inputs its range; but do it anyway... ;-))
out.pmaxima = length(dmaximai)/floor(lami/2);
out.modeperiodmax = mode(dmaximai);
out.pmodeperiodmax = sum(dmaximai == mode(dmaximai))/length(dmaximai);

% Same for local minima
% Look for periodicities in local maxima
minimai = find(dami(1:end-1) < 0 & dami(2:end) > 0) + 1;
dminimai = diff(minimai);
% is there a big peak in dmaxima?
% (no need to normalize since a given method inputs its range; but do it anyway... ;-))
out.pminima = length(dminimai)/floor(lami/2);
out.modeperiodmin = mode(dminimai);
out.pmodeperiodmin = sum(dminimai == mode(dminimai))/length(dminimai);

% Number of crossings at mean/median level, percentiles:
out.pcrossmean = sum(BF_sgnchange(ami-mean(ami)))/(lami-1);
out.pcrossmedian = sum(BF_sgnchange(ami-median(ami)))/(lami-1);
out.pcrossq10 = sum(BF_sgnchange(ami-quantile(ami,0.1)))/(lami-1);
out.pcrossq90 = sum(BF_sgnchange(ami-quantile(ami,0.9)))/(lami-1);

% ac1
out.amiac1 = CO_AutoCorr(ami,1,'Fourier');

end
