% ------------------------------------------------------------------------------
% NL_BoxCorrDim
% ------------------------------------------------------------------------------
% 
% References TSTOOL code, corrdim, to estimate the correlation dimension of a
% time-delay embedded time series using a box-counting approach.
% 
% TSTOOL: http://www.physik3.gwdg.de/tstool/
% 
%---INPUTS:
% y, column vector of time series data
% 
% nbins, maximum number of partitions per axis
% 
% embedparams [opt], embedding parameters as {tau,m} in 2-entry cell for a
%                   time-delay, tau, and embedding dimension, m. As inputs to BF_embed.
% 
%---OUTPUTS: Simple summaries of the outputs from corrdim.
% 
%---HISTORY:
% Ben Fulcher, November 2009
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = NL_BoxCorrDim(y,nbins,embedparams)

doPlot = 0; % plot outputs to a figure

% ------------------------------------------------------------------------------
%% Check inputs, preliminaries
% ------------------------------------------------------------------------------
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

% ------------------------------------------------------------------------------
%% Embed the signal
% ------------------------------------------------------------------------------
% convert to embedded signal object for TSTOOL
s = BF_embed(y,embedparams{1},embedparams{2},1);

if ~strcmp(class(s),'signal') && isnan(s); % embedding failed
    error('Time-series embedding to signal class for TSTOOL failed')
end

% ------------------------------------------------------------------------------
%% Run
% ------------------------------------------------------------------------------
rs = data(corrdim(s,nbins));

% Contains ldr as rows for embedding dimensions 1:m as columns;
if doPlot
    figure('color','w'); box('on');
    plot(rs,'k');
end

% ------------------------------------------------------------------------------
%% Output Statistics
% ------------------------------------------------------------------------------
% These statistics are just from my intuition

m = size(rs,2); % number of embedding dimensions
ldr = size(rs,1); % not completely clear from TSTOOL what ldr represents (= 17)

for i = 2:m
    out.(sprintf('meand%u',i)) = mean(rs(:,i));
    out.(sprintf('mediand%u',i)) = median(rs(:,i));
    out.(sprintf('mind%u',i)) = min(rs(:,i));
end

for i = 2:ldr
    out.(sprintf('meanr%u',i)) = mean(rs(i,:));
    out.(sprintf('medianr%u',i)) = median(rs(i,:));
    out.(sprintf('minr%u',i)) = min(rs(i,:));
    out.(sprintf('meanchr%u',i)) = mean(diff(rs(i,:)));
end

out.stdmean = std(mean(rs));
out.stdmedian = std(median(rs));

rsstretch = rs(:);
out.medianstretch = median(rsstretch);
out.minstretch = min(rsstretch);
out.iqrstretch = iqr(rsstretch);

end