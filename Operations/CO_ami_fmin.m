function out = CO_ami_fmin(y,meth,nbins)
% Analyzes patterns in the automutual information of a time series over a range of time lags up to the first minimum
% Automutual informations are calculated using a histogram method using CO_ami_benhist
% If nbins = [ . . . ] is a vector -- varies over this range of 
% Ben Fulcher September 2009

%% Set defaults:
% Default number of bins
if nargin < 3,
	nbins = 10;
end

N = length(y);

% Range of time lags, tau, to consider
%	(although loop usually broken before this maximum)
taur = 0:1:round(N/2);
ntaur = length(taur);

% Range of bin numbers to consider
nbinr = length(nbins);
amimins = zeros(nbinr,1);

% Calculate automutual information
for i = 1:nbinr % vary over number of bins in histogram
    amis = zeros(ntaur,1);
    for j = 1:ntaur % vary over time lags, tau
        amis(j) = CO_ami_benhist(y,taur(j),meth,nbins(i));
        if j > 2 && (amis(j)-amis(j-1))*(amis(j-1)-amis(j-2)) < 0
            amimins(i) = taur(j-1);
            break
        end
    end
    if amimins(i) == 0
		amimins(i) = taur(end)
	end
end

% plot(amimins,'o-k');

% Things to look for in the variation
% 1) Basic statistics
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