function out=CO_embed2_dist(y,tau)
% Looks at distance distributions in the two-dimensional
% embedding plot, and at time delay tau
% y should be z-scored
% Ben Fulcher September 2009
% Ben Fulcher 19/3/2010 -- fixed error in choosing tau too large

if nargin < 2 || isempty(tau)
    tau = 'tau';
end

if strcmp(tau,'tau'),
    tau = CO_fzcac(y);
    if tau>length(y)/10
        tau = floor(length(y)/10);
    end
end

if size(y,2)>size(y,1);
    y=y';
end

m=[y(1:end-tau) y(1+tau:end)];
% m2=y(1+tau:end);
N = size(m,1);

% plot(m(:,1),m(:,2),'.');
% input('this is your space!')

% distribution of euclidean distances between successive points in
% this space

% d = sum(abs(diff(m(:,1)).^2 + diff(m(:,1)).^2),2); % was like this until 9/7/2011 :-(
d = sum(abs(diff(m(:,1)).^2 + diff(m(:,2)).^2),2); % was like this from 9/7/2011 :-)

out.d_ac1 = CO_autocorr(d,1);
out.d_ac2 = CO_autocorr(d,2);
out.d_ac3 = CO_autocorr(d,3);
out.d_median = median(d);
out.d_max = max(d);
out.d_iqr = iqr(d);
out.d_min = min(d);
out.d_mean = mean(d);
out.d_std = std(d);
out.d_cv = mean(d)/std(d);


% distribution often fits exponential quite well
% fit to all values (often some extreme outliers, but ok)
% histogram with fixed bins maybe biased, but we'll see...
[n x] = hist(d,20);
n = n/(sum(n)*(x(2)-x(1))); % normalize
l = expfit(d);
nlogL = explike(l,d);
expf = exppdf(x,l);
out.d_expfit_l = l;
out.d_expfit_nlogL = nlogL;
 % sum of abs differences between exp fit and observed:
out.d_expfit_sumdiff = sum(abs(n - expf));

% keyboard


end