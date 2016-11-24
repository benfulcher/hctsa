function out = NL_TSTL_ReturnTime(y,NNR,maxT,past,Nref,embedParams)
% NL_TSTL_ReturnTime    Analysis of the histogram of return times.
%
% Return times are the time taken for the time series to return to a similar
% location in phase space for a given reference point
%
% Strong peaks in the histogram are indicative of periodicities in the data.
%
%---INPUTS:
%
% y, scalar time series as a column vector
% NNR, number of nearest neighbours
% maxT, maximum return time to consider
% past, Theiler window
% Nref, number of reference indicies
% embedParams, to feed into BF_embed
%
%---OUTPUTS: include basic measures from the histogram, including the occurrence of
% peaks, spread, proportion of zeros, and the distributional entropy.
%
% Uses the code, return_time, from TSTOOL.
% TSTOOL: http://www.physik3.gwdg.de/tstool/

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
%% Check Inputs
% ------------------------------------------------------------------------------
N = length(y); % length of the input time series

% Number of nearest neighbours, NNR
if nargin < 2 || isempty(NNR)
    NNR = 5;
end
if (NNR > 0) && (NNR < 1) % specify a proportion of time series length
    NNR = floor(NNR*N); if NNR == 0, NNR = 1; end
end

% Maximum return time, maxT
if nargin < 3 || isempty(maxT)
    maxT = 0.1;
end
if (maxT > 0) && (maxT <= 1) % specify a proportion
    maxT = floor(N*maxT);
    if maxT == 0, maxT = 1; end
end

% Theiler window, past
if nargin < 4 || isempty(past)
    past = 10;
end
if (past > 0) && (past < 1) % specify a proportion
    past = floor(N*past);
    if past == 0, past = 1; end % round up from 0
end

% Number of reference points
if nargin < 5 || isempty(Nref)
    Nref = -1; % use all available points
end

% embed parameters
if nargin < 6 || isempty(embedParams)
    embedParams = {'ac','fnnmar'};
    fprintf(1,'Using default embedding using autocorrelation and cao\n');
end

doPlot = 0; % plot outputs to figures

% ------------------------------------------------------------------------------
%% Embed the signal
% ------------------------------------------------------------------------------
s = BF_embed(y,embedParams{1},embedParams{2},1);

if ~isa(s,'signal') && isnan(s); % embedding failed
    fprintf(1,'Embedding failed\n');
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Run the code
% ------------------------------------------------------------------------------
try
    rs = return_time(s, NNR, maxT, past, Nref);
catch emsg
    if strcmp(emsg.message,'Index exceeds matrix dimensions.')
        fprintf(1,'Error evaluating return_time\n');
        out = NaN; return
    else
        error(emsg.message);
    end
end

Trett = data(rs);
NN = length(Trett);

% ------------------------------------------------------------------------------
%% Quantify structure in output
% ------------------------------------------------------------------------------
out.max = max(Trett);
out.std = std(Trett);
out.pzeros = sum(Trett == 0)/NN;
out.pg05 = sum(Trett>max(Trett)*0.5)/NN;
out.iqr = iqr(Trett);

% recurrent peaks:
icross05 = find((Trett(1:end-1)-0.5*max(Trett)).*(Trett(2:end)-0.5*max(Trett)) < 0);
if ~isempty(icross05) && length(icross05) > 2
    difficross05 = diff(icross05);
    difficross05 = difficross05(difficross05 > 0.4*max(difficross05)); % remove small entries, crossing peaks

    out.meanpeaksep = mean(difficross05)/NN;
    out.maxpeaksep = max(difficross05)/NN;
    out.minpeaksep = min(difficross05)/NN;
    out.rangepeaksep = range(difficross05)/NN;
    out.stdpeaksep = std(difficross05)/sqrt(NN);
else
    out.meanpeaksep = NaN;
    out.maxpeaksep = NaN;
    out.minpeaksep = NaN;
    out.rangepeaksep = NaN;
    out.stdpeaksep = NaN;
end

out.statrtys = std(Trett(1:floor(end/2)))/std(Trett(floor(end/2)+1:end));
out.statrtym = mean(Trett(1:floor(end/2)))/mean(Trett(floor(end/2)+1:end));

out.hhist = -sum(Trett(Trett>0).*log(Trett(Trett>0)));

% ------------------------------------------------------------------------------
%% Coarse-grain to 20 bins
% ------------------------------------------------------------------------------
numBins = 20;
cglav = zeros(numBins,1);
inds = round(linspace(0,NN,numBins+1));
for i = 1:numBins
    cglav(i) = sum(Trett(inds(i)+1:inds(i+1)));
end
if doPlot
    figure('color','w'); box('on'); plot(cglav,'k')
end
out.hcgdist = -sum(cglav(cglav > 0).*log(cglav(cglav > 0)));
out.rangecgdist = range(cglav);
out.pzeroscgdist = sum(cglav == 0)/numBins;

% ------------------------------------------------------------------------------
%% Get distribution of distribution of return times
% ------------------------------------------------------------------------------
[nhist, binEdges] = histcounts(Trett,'BinMethod','sqrt','Normalization','probability');
if doPlot
    binCenters = mean([binEdges(1:end-1); binEdges(2:end)]);
    figure('color','w'); plot(binCenters,nhist,'o-k')
end
out.maxhisthist = max(nhist);
out.phisthistmin = nhist(1); % this is the same as maxhisthist
out.hhisthist = -sum(nhist(nhist > 0).*log(nhist(nhist > 0)));

end
