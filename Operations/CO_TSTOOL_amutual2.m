% CO_TSTL_amutual2
% 
% Uses amutual2 code from TSTOOL up to a maximum time-delay maxtau.
% 
% TSTOOL: http://www.physik3.gwdg.de/tstool/
% 
% y, column vector of time series data
% maxtau, maximal lag
% 
% Statistics are returned on the output of amutual2 over this range, as for
% CO_TSTL_amutual.
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

function out = CO_TSTL_amutual2(y,maxtau)
% Ben Fulcher, 2009

%% Preliminaries
N = length(y); % length of time series
s = signal(y); % convert to signal object for TSTOOL

%% Check Inputs
maxtau0 = maxtau;
if nargin < 2 || isempty(maxtau)
    maxtau = ceil(N/4);
else
    maxtau = min(maxtau,ceil(N/2)); % if given a number, don't go above N/2
end

%% Run
ami = data(amutual2(s,maxtau));

%% Give output
% (c.f., CO_TSTL_amutual -- similar routine here)

% output the raw values
for i = 1:maxtau0
    if i <= maxtau
        eval(sprintf('out.ami%u = ami(%u);',i,i));
    else
        eval(sprintf('out.ami%u = NaN;',i));
    end
end

% output statistics
lami = length(ami);

% mean mutual information over this lag range
out.mami = mean(ami);
out.stdami = std(ami);


% first miniimum of mutual information across range
dami = diff(ami);
extremai = find(dami(1:end-1).*dami(2:end) < 0);
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

% Same for local minima
% Look for periodicities in local maxima
minimai = find(dami(1:end-1) < 0 & dami(2:end) > 0) + 1;
dminimai = diff(minimai);
% is there a big peak in dmaxima?
 % (no need to normalize since a given method inputs its range; but do it anyway... ;-))
out.pminima = length(dminimai)/floor(lami/2);
out.modeperiodmin = mode(dminimai);
out.pmodeperiodmin = sum(dminimai == mode(dminimai))/length(dminimai);

% number of crossings at mean/median level, percentiles
out.pcrossmean = sum(BF_sgnchange(ami-mean(ami)))/(lami-1);
out.pcrossmedian = sum(BF_sgnchange(ami-median(ami)))/(lami-1);
out.pcrossq10 = sum(BF_sgnchange(ami-quantile(ami,0.1)))/(lami-1);
out.pcrossq90 = sum(BF_sgnchange(ami-quantile(ami,0.9)))/(lami-1);

% ac1
out.amiac1 = CO_AutoCorr(ami,1);

end