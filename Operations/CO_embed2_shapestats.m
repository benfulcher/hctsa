function out = CO_embed2_shapestats(y,tau,shape,r)
% Takes a shape around with the points and counts the points inside it
% This becomes a time series that can be analyzed... (there's a temporal verison as well, that
% does this for the time series plotted as a function of time...): CO_t_shape_translate
% The input time series, y, should be a z-scored column vector
% Implements using circles of radius r
% Ben Fulcher September 2009
% Ben Fulcher 19/3/2010 -- fixed error in choosing tau too high for highly
%                           autocorrelated signals using option 'tau'

doplot = 0; % plot results for debugging

%% Check inputs, set defaults:
if nargin < 2 || isempty(tau)
	fprintf(1,'Using tau as first zero crossing ' ...
			'of the autocorrelation function.\n')
    tau = 'tau';
end
if nargin < 3 || isempty(shape)
	shape = 'circle';
	fprintf(1,'Using a circle.\n');
end
if nargin < 4 || isempty(r)
	
end

% Can set time lag equal to first zero crossing of the autocorrelation function with the 'tau' input
if strcmp(tau,'tau'),
    tau = CO_fzcac(y);
    if tau > length(y)/10
        tau = floor(length(y)/10);
    end
end

% Ensure y is a column vector:
if size(y,2) > size(y,1);
	y = y';
end

%% Create the recurrence space, populated by points m
m = [y(1:end-tau), y(1+tau:end)];
N = length(m);

if doplot % plot the recurrence space:
	plot(m(:,1),m(:,2),'.');
end

%% Start the analysis

counts = zeros(N,1); % stores the counts for points within circle

switch shape
	case 'circle'
		% Puts a circle around each point in the embedding space in turn
		% counts how many points are inside this shape, looks at the time series thus formed

		for i = 1:N % across all points in the time series
			
			m_c = m - ones(N,1)*m(i,:); % points wrt current point i
			m_c_d = sum(m_c.^2,2); % Euclidean distances from point i
			
		    counts(i) = sum(m_c_d <= r^2); % number of points enclosed in a circle of radius r
		end
    otherwise
        error('Unknown shape ''%s''', shape)
end
counts = counts - 1; % ignore self-counts

% No meaningful output if never got a count with any other point!
% (radius, r, is probably too small)
if all(counts == 0)
    fprintf(1,'No counts detected!\n');
    out = NaN; return
end

% Return basic statistics on the counts
out.ac1 = CO_autocorr(counts,1);
out.ac2 = CO_autocorr(counts,2);
out.ac3 = CO_autocorr(counts,3);
out.tau = CO_fzcac(counts);
out.max = max(counts);
out.std = std(counts);
out.median = median(counts);
out.iqr = iqr(counts);
out.iqronrange = out.iqr/range(counts);

% Distribution
x = (0:max(counts));
n = hist(counts,x); n = n/sum(n);
[out.mode_val, mix] = max(n);
out.mode = x(mix);

% Poisson fit to distribution
% (note that this is actually rarely a good fit...)
l = poissfit(counts);
poiss_n = poisspdf(x,l);
out.poissfit_l = l;
% goodness of fit:
out.poissfit_absdiff = sum(abs(poiss_n-n));
out.poissfit_sqdiff = sum((poiss_n-n).^2);

% Entropy in 10-bin histogram
[n, x] = hist(counts,10); n = n/(sum(n)*(x(2)-x(1)));
out.hist10_ent = sum(n(n>0).*log(n(n > 0)));

% plot(x,poisspdf(x,l),'g'); hold on;
% plot(x,n,'k'); hold off

% Stationarity measure for fifths of the time series
afifth = floor(N/5);
buffer_m = zeros(afifth,5); % stores a fifth of the time series (embedding vector) in each entry
for i = 1:5
	buffer_m(:,i) = counts(afifth*(i-1)+1:afifth*i);
end
out.statav5_m = std(mean(buffer_m))/std(counts);
out.statav5_s = std(std(buffer_m))/std(counts);

end