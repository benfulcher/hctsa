function out = CO_ami_dynbenhist(y,meth,nbins)
% Analyzes the automutual information of a time series using histograms computed using CO_ami_benhist
% The input time series, y, should be z-scored
% The automutual information is computed using a specified number of bins, nbins
% Ben Fulcher September 2009

%% Set defaults:
% Default number of bins
if nargin < 3,
	nbins = 10;
end

maxtau = 50; % maximum time lag to compute up to

%% Do the calculation:
taur = 0:1:maxtau; % define the set of time lags
nr = length(taur); % the total number of lags to compute the AMI at
amis = zeros(nr,1); % the automutual information as a function of time lag

for i = 1:nr;
    amis(i) = CO_ami_benhist(y,taur(i),meth,nbins);
end

% plot(taur,amis,'.-k');


end