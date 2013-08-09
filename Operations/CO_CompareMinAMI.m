% CO_CompareMinAMI
% 
% Finds the first minimum of the automutual information by various different
% estimation methods, and sees how this varies over different coarse-grainings
% of the time series.
% 
% The function returns a set of statistics on the set of first minimums of the
% automutual information function obtained over a range of the number of bins
% used in the histogram estimation, when specifying 'nbins' as a vector
% 
% INPUTS:
% y, the input time series
% 
% meth, the method for estimating mutual information (input to CO_HistogramAMI)
% 
% nbins, the number of bins for the AMI estimation to compare over (can be a
%           scalar or vector)
% 
% Outputs include the minimum, maximum, range, number of unique values, and the
% position and periodicity of peaks in the set of automutual information
% minimums.
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

function out = CO_CompareMinAMI(y,meth,nbins)
% Ben Fulcher, September 2009

%% Set defaults:
% Default number of bins
if nargin < 3,
	nbins = 10;
end

doplot = 0; % plot outputs to figure
N = length(y); % time-series length

% Range of time lags, tau, to consider
%	(although loop usually broken before this maximum)
taur = (0:1:round(N/2));
ntaur = length(taur);

% Range of bin numbers to consider
nbinr = length(nbins);
amimins = zeros(nbinr,1);

% Calculate automutual information
for i = 1:nbinr % vary over number of bins in histogram
    amis = zeros(ntaur,1);
    for j = 1:ntaur % vary over time lags, tau
        amis(j) = CO_HistogramAMI(y,taur(j),meth,nbins(i));
        if (j > 2) && ((amis(j)-amis(j-1))*(amis(j-1)-amis(j-2)) < 0)
            amimins(i) = taur(j-1);
            break
        end
    end
    if amimins(i) == 0
		amimins(i) = taur(end)
	end
end

if doplot
    figure('color','w');
    plot(amimins,'o-k');
end

% Things to look for in the variation
% Basic statistics
out.min = min(amimins);
out.max = max(amimins);
out.range = range(amimins);
out.median = median(amimins);
out.mean = mean(amimins);
out.std = std(amimins);
out.iqr = iqr(amimins);

% Unique values, mode
out.nunique = length(unique(amimins));
[out.mode, out.modef] = mode(amimins);
out.modef = out.modef/nbinr;
% hist = zeros(length(u),1);
% for i=1:length(u)
%     hist(i) = sum(n == u(i));
% end
% out.mode = u(find(hist == max(hist),1,'first'));

% Converged value?
out.conv4 = mean(amimins(end-4:end));

% Look for peaks
% local maxima above 1*std from mean
% inspired by curious result of periodic maxima for periodic signal with
% bin size... ('quantiles', [2:80])
loc_extr = intersect(find(diff(amimins(1:end-1)) > 0), BF_sgnchange(diff(amimins(1:end-1)),1)) + 1;
big_loc_extr = intersect(find(amimins > out.mean+out.std),loc_extr);
out.nlocmax = length(big_loc_extr);
if out.nlocmax > 2
    out.maxp = std(diff(big_loc_extr));
else
    out.maxp = NaN;
end


end