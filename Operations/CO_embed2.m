function out = CO_embed2(y,tau)
% Analyzes angular distributions and other such statistics in a two-dimensional recurrence plot at a time delay tau
% The input time series, y, should be a z-scored column vector
% Ben Fulcher September 2009
% Ben Fulcher 19/3/2010 -- corrected error in choosing tau too large with
%                           'tau'

doplot = 0; % can set to 1 to plot some outputs

%% Set defaults
if nargin < 2 || isempty(tau)
    tau = 'tau';
end

% Set tau to the first zero-crossing of the autocorrelation function, with the 'tau' input
if strcmp(tau,'tau'),
    tau = CO_fzcac(y);
    if tau > length(y)/10
        tau = floor(length(y)/10);
    end
end

% Ensure that y is a column vector
if size(y,2) > size(y,1);
  y = y';
end

% Construct the two-dimensional recurrence space
m = [y(1:end-tau), y(1+tau:end)];
N = size(m,1); % number of points in the recurrence space

if doplot
  figure('color','w'); plot(m(:,1),m(:,2),'.');
end

% 1) Distribution of angles time series; angles between successive points in
% 	 this space

theta = diff(m(:,2))./diff(m(:,1));
theta = atan(theta); % measured as deviation from the horizontal


if doplot, ksdensity(theta); end % can plot distribution of angles
out.theta_ac1 = CO_autocorr(theta,1);
out.theta_ac2 = CO_autocorr(theta,2);
out.theta_ac3 = CO_autocorr(theta,3);

out.theta_mean = mean(theta);
out.theta_std = std(theta);

x = linspace(-pi/2,pi/2,11); % 10 bins in the histogram
n = histc(theta,x); n(end-1)=n(end-1)+n(end); n=n(1:end-1); n=n/sum(n);
out.hist10std = std(n);
out.histent = -sum(n(n>0).*log(n(n>0)));

% Stationarity in fifths of the time series
% Use histograms with 4 bins
x = linspace(-pi/2,pi/2,5); % 4 bins
afifth = floor((N-1)/5); % -1 because angles are correlations *between* points
n = zeros(length(x),5);
for i = 1:5
	n(:,i) = histc(theta(afifth*(i-1)+1:afifth*i),x);
end
n = n/afifth;
n(4,:) = n(4,:) + n(5,:); n(5,:) = [];

% Output the standard deviation in each bin
out.stdb1 = std(n(:,1));
out.stdb2 = std(n(:,2));
out.stdb3 = std(n(:,3));
out.stdb4 = std(n(:,4));


% Points in the space
% Stationarity of points in the space (do they move around in the space)

% (1) in terms of distance from origin
afifth = floor(N/5);
buffer_m = cell(5,1); % stores a fifth of the time series (embedding vector) in each entry
for i = 1:5
	buffer_m{i} = m(afifth*(i-1)+1:afifth*i,:);
end

% Mean euclidean distance in each segment
eucdm = cellfun(@(x)mean(sqrt(x(:,1).^2 + x(:,2).^2)),buffer_m);
out.eucdm1 = eucdm(1); out.eucdm2 = eucdm(2); out.eucdm3 = eucdm(3);
out.eucdm4 = eucdm(4); out.eucdm5 = eucdm(5);
out.std_eucdm = std(eucdm); out.mean_eucdm = mean(eucdm);

% Standard deviation of Euclidean distances in each segment
eucds = cellfun(@(x)std(sqrt(x(:,1).^2 + x(:,2).^2)),buffer_m);
out.eucds1 = eucds(1); out.eucds2 = eucds(2); out.eucds3 = eucds(3);
out.eucds4 = eucds(4); out.eucds5 = eucds(5);
out.std_eucds = std(eucds); out.mean_eucds = mean(eucds);

% Maximum volume in each segment
% (defined as area of rectangle of max span in each direction)
maxspanx = cellfun(@(x)range(x(:,1)),buffer_m);
maxspany = cellfun(@(x)range(x(:,2)),buffer_m);
spanareas = maxspanx.*maxspany;
out.stdspana = std(spanareas);
out.meanspana = mean(spanareas);

% Outliers
% area of max span of all points; versus area of max span of 50% of points closest to origin
d = sqrt(m(:,1).^2 + m(:,2).^2);
[d_sort, ix] = sort(d,'ascend');

out.areas_all = range(m(:,1))*range(m(:,2));
r50 = ix(1:round(end/2)); % 50% of point closest to origin
out.areas_50 = range(m(r50,1))*range(m(r50,2));
out.arearat = out.areas_50 / out.areas_all;

end