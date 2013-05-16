function out = CO_embed2_shapestats(y,tau,shape,r)
% Takes a shape around with the points and counts the points inside it
% This becomes a time series that can be analyzed... (there's a temporal verison as well, that
% does this for the time series plotted as a function of time...): CO_t_shape_translate
% y should be z-scored and a column vector
% Implements with a cirle or radius r
% Ben Fulcher September 2009
% Ben Fulcher 19/3/2010 -- fixed error in choosing tau too high for highly
%                           autocorrelated signals using option 'tau'

%% Preliminaries
if nargin<2 || isempty(tau)
    tau = 'tau';
end

if strcmp(tau,'tau'),
    tau = CO_fzcac(y);
    if tau > length(y)/10
        tau = floor(length(y)/10);
    end
end

if size(y,2) > size(y,1); y = y'; end

%% Create the embedding space
m = [y(1:end-tau) y(1+tau:end)];
N = length(m);
% plot(m(:,1),m(:,2),'.');
% input('this is your space!')


%% Do your thing

switch shape
	case 'circle'
		% puts a circle around each point in the embedding space in turn
		% counts how many points are inside this shape, looks at the time series thus formed

		silj = zeros(N,1); % stores the statistics

		for i = 1:N % across all points in the time series
			
			m_c = m - ones(N,1)*m(i,:); % points wrt current point i
			m_c_d = sum(m_c.^2,2); % euclidean distances from point i
			
		    silj(i) = length(find(m_c_d <= r^2)); % number of points enclosed in circle
		
			% end
		end
end
silj = silj-1; % ignore self-counts

if all(silj==0)
    out = NaN; return
end

% return basic statistics on this silj
out.ac1 = CO_autocorr(silj,1);
out.ac2 = CO_autocorr(silj,2);
out.ac3 = CO_autocorr(silj,3);
out.tau = CO_fzcac(silj);
out.max = max(silj);
out.std = std(silj);
out.median = median(silj);
out.iqr = iqr(silj);
out.iqronrange = out.iqr/range(silj);

% distribution
x = 0:max(silj);
n = hist(silj,x); n=n/sum(n);
[out.mode_val mix] = max(n);
out.mode = x(mix);

% poisson fit to distribution
% (note that this is actually rarely a good fit...)
l = poissfit(silj);
poiss_n = poisspdf(x,l);
out.poissfit_l = l;
% gof
out.poissfit_absdiff = sum(abs(poiss_n-n));
out.poissfit_sqdiff = sum((poiss_n-n).^2);

% entropy in 10-bin histogram
[n x] = hist(silj,10); n=n/(sum(n)*(x(2)-x(1)));
out.hist10_ent = sum(n(n>0).*log(n(n>0)));

% plot(x,poisspdf(x,l),'g'); hold on;
% plot(x,n,'k'); hold off

% Stationarity 5
afifth=floor(N/5);
buffer_m=zeros(afifth,5); % stores a fifth of the time series (embedding vector) in each entry
for i=1:5
	buffer_m(:,i) = silj(afifth*(i-1)+1:afifth*i);
end
out.statav5_m = std(mean(buffer_m))/std(silj);
out.statav5_s = std(std(buffer_m))/std(silj);

% keyboard


end