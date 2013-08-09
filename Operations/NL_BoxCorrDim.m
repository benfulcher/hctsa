% NL_BoxCorrDim
% 
% References TSTOOL code, corrdim, to estimate the correlation dimension of a
% time-delay embedded time series using a box-counting approach.
% 
% TSTOOL: http://www.physik3.gwdg.de/tstool/
% 
% INPUTS:
% y, column vector of time series data
% 
% nbins, maximum number of partitions per axis
% 
% embedparams [opt], embedding parameters as {tau,m} in 2-entry cell for a
%                   time-delay, tau, and embedding dimension, m. As inputs to BF_embed.
% 
% Output statistics are simple summaries of the outputs from this algorithm.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
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
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = NL_BoxCorrDim(y,nbins,embedparams)
% Ben Fulcher, November 2009

doplot = 0; % plot outputs to a figure

%% Preliminaries
N = length(y); % length of time series

% (1) Maxmum number of partitions per axis, nbins
if nargin < 2 || isempty(nbins)
    nbins = 100; % default number of bins per axis is 100
end

% (2) Set embedding parameters to defaults
if nargin < 3 || isempty(embedparams)
    embedparams = {'ac','cao'};
else
    if length(embedparams) ~= 2
        error('Embedding parameters should be formatted like {tau,m}')
    end
end

%% Embed the signal
% convert to embedded signal object for TSTOOL
s = BF_embed(y,embedparams{1},embedparams{2},1);

if ~strcmp(class(s),'signal') && isnan(s); % embedding failed
    error('Time-series embedding to signal class for TSTOOL failed')
end

%% Run
rs = data(corrdim(s,nbins));
% Contains ldr as rows for embedding dimensions 1:m as columns;
if doplot
    figure('color','w'); box('on');
    plot(rs,'k');
end

%% Output Statistics
% Note: these aren't particularly well motivated.
m = size(rs,2); % number of embedding dimensions
ldr = size(rs,1); % I don't really know what this means; = 17
for i = 2:m
    meani = mean(rs(:,i));
    eval(sprintf('out.meand%u = meani;',i))
    mediani = median(rs(:,i));
    eval(sprintf('out.mediand%u = mediani;',i))
    mini = min(rs(:,i));
    eval(sprintf('out.mind%u = mini;',i))
end

for i = 2:ldr
    meani = mean(rs(i,:));
    eval(sprintf('out.meanr%u = meani;',i))
    mediani = median(rs(i,:));
    eval(sprintf('out.medianr%u = mediani;',i))
    mini = min(rs(i,:));
    eval(sprintf('out.minr%u = mini;',i))
    meanchi = mean(diff(rs(i,:)));
    eval(sprintf('out.meanchr%u = meanchi;',i))
end

out.stdmean = std(mean(rs));
out.stdmedian = std(median(rs));

rsstretch = rs(:);
out.medianstretch = median(rsstretch);
out.minstretch = min(rsstretch);
out.iqrstretch = iqr(rsstretch);


end