function out=CO_embed2(y,tau)
% Looks at angular distributions and other such statistics in the two-dimensional
% recurrencey sort of plot thing, and at time delay tau
% y should be z-scored
% Ben Fulcher September 2009
% Ben Fulcher 19/3/2010 -- corrected error in choosing tau too large with
%                           'tau'

if nargin<2 || isempty(tau)
    tau = 'tau';
end

if strcmp(tau,'tau'),
    tau = CO_fzcac(y);
    if tau > length(y)/10
        tau = floor(length(y)/10);
    end
end

if size(y,2)>size(y,1); y=y'; end

m = [y(1:end-tau) y(1+tau:end)];
% m2=y(1+tau:end);
N = size(m,1);

% plot(m(:,1),m(:,2),'.');
% input('this is your space!')

% 1) distribution of angles time series; angles between successive points in
% 	 this space

% theta=zeros(N,1);
theta = diff(m(:,2))./diff(m(:,1));
theta = atan(theta); % measured as deviation from the horizontal


% ksdensity(theta)
out.theta_ac1 = CO_autocorr(theta,1);
out.theta_ac2 = CO_autocorr(theta,2);
out.theta_ac3 = CO_autocorr(theta,3);

out.theta_mean = mean(theta);
out.theta_std = std(theta);

x = linspace(-pi/2,pi/2,11); % ten bins
n = histc(theta,x); n(end-1)=n(end-1)+n(end); n=n(1:end-1); n=n/sum(n);
out.hist10std = std(n);
out.histent=-sum(n(n>0).*log(n(n>0)));

% stationarity in 5; histograms with 4 bins
x=linspace(-pi/2,pi/2,5); % 4 bins
afifth=floor((N-1)/5); % -1 because angles are correlations *between* points
n=zeros(length(x),5);
for i = 1:5
	n(:,i) = histc(theta(afifth*(i-1)+1:afifth*i),x);
end
n = n/afifth;
n(4,:) = n(4,:) + n(5,:); n(5,:) = [];

% standard deviation in each bin
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

% mean euclidean distance in each segment
eucdm = cellfun(@(x)mean(sqrt(x(:,1).^2+x(:,2).^2)),buffer_m);
out.eucdm1 = eucdm(1); out.eucdm2 = eucdm(2); out.eucdm3 = eucdm(3);
out.eucdm4 = eucdm(4); out.eucdm5 = eucdm(5);
out.std_eucdm = std(eucdm); out.mean_eucdm = mean(eucdm);

% std of euclidean distances in each segment
eucds = cellfun(@(x)std(sqrt(x(:,1).^2+x(:,2).^2)),buffer_m);
out.eucds1 = eucds(1); out.eucds2 = eucds(2); out.eucds3 = eucds(3);
out.eucds4 = eucds(4); out.eucds5 = eucds(5);
out.std_eucds = std(eucds); out.mean_eucds = mean(eucds);

% max volume in each segment
% (defined as area of rectangle of max span in each direction)
maxspanx = cellfun(@(x)range(x(:,1)),buffer_m);
maxspany = cellfun(@(x)range(x(:,2)),buffer_m);
spanareas = maxspanx.*maxspany;
out.stdspana = std(spanareas);
out.meanspana = mean(spanareas);

% Outliers
% area of max span of all points; versus area of max span of 50% of points closest to origin
d = sqrt(m(:,1).^2+m(:,2).^2);
[d_sort ix] = sort(d,'ascend');

out.areas_all = range(m(:,1))*range(m(:,2));
r50 = ix(1:round(end/2)); % 50% of point closest to origin
out.areas_50 = range(m(r50,1))*range(m(r50,2));
out.arearat = out.areas_50 / out.areas_all;


end