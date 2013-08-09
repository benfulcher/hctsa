% CO_TSTL_amutual
% 
% Uses amutual code from TSTOOL, which uses a
% histogram method with n bins to estimate the mutual information of a
% time series across a range of time-delays, tau.
% 
% TSTOOL: http://www.physik3.gwdg.de/tstool/
%
% INPUTS:
% 
% y, the time series
% 
% maxtau, the maximum lag for which to calculate the auto mutual information
% 
% nbins, the number of bins for histogram calculation
% 
% A number of statistics of the function over the range of tau are returned,
% including the mean mutual information, its standard deviation, first minimum,
% proportion of extrema, and measures of periodicity in the positions of local
% maxima.
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

function out = CO_TSTL_amutual(y,maxtau,nbins)
% Ben Fulcher, October 2009

doplot = 0; % toggle plotting of outputs

%% Preliminaries
N = length(y); % length of time series
s = signal(y); % convert to signal object for TSTOOL

%% Existence checks
if ~exist('amutual')
    error('''amutual'' not found -- ensure the TSTOOL package is installed correctly??\n');
end

%% Check Inputs
if nargin < 2 || isempty(maxtau)
    maxtau = ceil(N/4);
end

if nargin < 3 || isempty(nbins)
    nbins = round(sqrt(N/10)); % this is an arbitrary choice (!!) ;-)
end

%% Run
ami = data(amutual(s,maxtau,nbins));
lami = length(ami);

% Plot results
if doplot
    figure('color','w'); box('on');
    plot(ami,'-ok');
end

% Change automutual information vector to a structure for output
for i = 1:maxtau+1
    eval(sprintf('out.ami%u = ami(%u);',i,i));
end

% mean mutual information over this lag range
out.mami = mean(ami);
out.stdami = std(ami);

% first miniimum of mutual information across range
dami = diff(ami);
extremai = find(dami(1:end-1).*dami(2:end)<0);
out.pextrema = length(extremai)/(lami-1);
if isempty(extremai)
   out.fmmi = lami; % actually represents lag, because indexes don't but diff delays by 1
else
    out.fmmi = extremai(1);
end

% Look for periodicities in local maxima
maximai = find(dami(1:end-1) > 0 & dami(2:end) < 0) + 1;
dmaximai = diff(maximai);
% is there a big peak in dmaxima?
% (no need to normalize since a given method inputs its range; but do it anyway... ;-))
out.pmaxima = length(dmaximai)/floor(lami/2);
out.modeperiodmax = mode(dmaximai);
out.pmodeperiodmax = sum(dmaximai == mode(dmaximai))/length(dmaximai);

if doplot
    hold on;
    plot(maximai,ami(maximai),'or');
    hold off
end


end