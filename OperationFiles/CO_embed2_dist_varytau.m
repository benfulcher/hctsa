function out=CO_embed2_dist_varytau(y,maxtau)
% looks at how properties in two-dimensional embedding space change as vary tau
% This isn't very informative!!! SCRAP!!
% Ben Fulcher September 2009

taur=1:maxtau;
ntaur = length(taur);

if size(y,2)>size(y,1); y=y'; end

stats_store = zeros(7,ntaur);

for i=1:ntaur

	tau=taur(i);
	
	m=[y(1:end-tau) y(1+tau:end)];
	d = sum(abs(diff(m(:,1)).^2 + diff(m(:,1)).^2),2);

	stats_store(1,i)=CO_autocorr(d,1);
	stats_store(2,i)=median(d);
	stats_store(3,i)=max(d);
	stats_store(4,i)=iqr(d);
	stats_store(5,i)=min(d);
	stats_store(6,i)=mean(d);
	stats_store(7,i)=std(d);

end

% plot(stats_store');
% keyboard

% out.ac1_thetaac1 = CO_autocorr(stats_store(1,:),1);
% out.ac1_thetaac2 = CO_autocorr(stats_store(2,:),1);
% out.ac1_thetaac3 = CO_autocorr(stats_store(3,:),1);
% out.mean_thetaac1 = mean(stats_store(1,:));
% out.max_thetaac1 = max(stats_store(1,:));
% out.min_thetaac1 = min(stats_store(1,:));
% out.mean_thetaac2 = mean(stats_store(2,:));
% out.max_thetaac2 = max(stats_store(2,:));
% out.min_thetaac2 = min(stats_store(2,:));
% out.mean_thetaac3 = mean(stats_store(3,:));
% out.max_thetaac3 = max(stats_store(3,:));
% out.min_thetaac3 = min(stats_store(3,:));
% out.meanrat_thetaac12 = out.mean_thetaac1/out.mean_thetaac2;
% out.diff_thetaac12 = sum(abs(stats_store(2,:)-stats_store(1,:)));

% 
% 
% function out=CO_embed2_dist(y,tau)
% % Looks at distance distributions in the two-dimensional
% % embedding plot, and at time delay tau
% % y should be z-scored
% % Ben Fulcher September 2009
% 
% if strcmp(tau,'tau'), tau=CO_fzcac(y); end
% if size(y,2)>size(y,1); y=y'; end
% 
% 
% % m2=y(1+tau:end);
% N=length(m);
% 
% plot(m(:,1),m(:,2),'.');
% input('this is your space!')
% 
% % distribution of euclidean distances between successive points in
% % this space
% 
% d = sum(abs(diff(m(:,1)).^2 + diff(m(:,1)).^2),2);
% 
% 
% 
% 
% % distribution often fits exponential quite well
% % fit to all values (often some extreme outliers, but ok)
% % histogram with fixed bins maybe biased, but we'll see...
% [n x] = hist(d,20);
% n=n/(sum(n)*(x(2)-x(1))); % normalize
% l=expfit(d);
% nlogL = explike(l,d);
% expf = exppdf(x,l);
% out.d_expfit_l = l;
% out.d_expfit_nlogL = nlogL;
%  % sum of abs differences between exp fit and observed:
% out.d_expfit_sumdiff = sum(abs(n - expf));
% 
% 
% keyboard
% 
% 
% end





end