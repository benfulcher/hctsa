% CO_ami_benhist
% 
% Calculates the automutual information using histograms, using a given approach
% to binning the data.
% 
% Uses hist2.m function (renamed hist2_NK.m here) by Nedialko Krouchev, obtained
% from Matlab Central,
% http://www.mathworks.com/matlabcentral/fileexchange/12346-hist2-for-the-people
% [[hist2 for the people by Nedialko Krouchev, 20 Sep 2006 (Updated 21 Sep 2006)]]
% 
% The automutual information is calculated using:
% (i) 'even': evenly-spaced bins through the range of the time series,
% (ii) 'std1', 'std2': bins that extend only up to a multiple of the standard deviation from the
% mean of the time series to exclude outliers, or
% (iii) 'quantiles': equiprobable bins chosen using quantiles.
% 
% A time-lag can be specified as tau
% 
% Some methods, meth, require the extra parameter nbins


function out = CO_ami_benhist(y,tau,meth,nbins)
% Ben Fulcher, September 2009


%% INPUTS
% Time-lag, tau
if nargin < 2 || isempty(tau)
    tau = 1;  % time-lag of 1
end
if strcmp(tau,'tau')
    tau = CO_fzcac(y);
    fprintf(1,'tau = %u set to fist zero-crossing of ACF\n',tau);
end

if nargin < 3 || isempty(meth)
    meth = 'even'; % default
end

if nargin < 4
    nbins = 10;
end

% 1) Form the time-delay vectors y1 and y2
y1 = y(1:end-tau);
y2 = y(1+tau:end);

% Number of options:
% remove outliers first?, number of bins, range of bins, bin sizes

%% Bins by standard deviation (=1)

% same for both -- assume same distribution (true for stationary processes,
% or small lags)
switch meth
    case 'even'
        b = linspace(min(y)-0.1,max(y)+0.1,nbins+1); % +0.1 to make sure all points included
        
    case 'std1' % std bins up to 1
        b = linspace(-1,1,nbins+1);
        if min(y) < -1; b = [min(y)-0.1, b]; end
        if max(y) > 1; b = [b, max(y)+0.1]; end
            
    case 'std2'
        b = linspace(-2,2,nbins+1);
        if min(y) < -2; b = [min(y)-0.1, b]; end
        if max(y) > 2; b = [b, max(y)+0.1]; end
            
    case 'quantiles' % use quantiles with ~equal number in each bin
        b = quantile(y,linspace(0,1,nbins+1));
        b(1) = b(1) - 0.1; b(end) = b(end) + 0.1;
        
    otherwise
        error('Unknown method ''%s''',meth)
end

nb = length(b) - 1; % number of bins (-1 since b defines edges)

% (1) Joint distribution of y1 and y2
pij = hist2_NK(y1,y2,b,b);
pij = pij(1:nb,1:nb); % joint
pij = pij/sum(sum(pij)); % joint
pi = sum(pij,1); % marginal
pj = sum(pij,2); % other marginal
% Old-fashioned method (should give same result):
% pi = histc(y1,b); pi = pi(1:nb); pi = pi/sum(pi); % marginal
% pj = histc(y2,b); pj= pj(1:nb); pj = pj/sum(pj); % other marginal

pii = ones(nb,1)*pi;
pjj = pj*ones(1,nb);

r = (pij > 0); % defining the range in this way, we set log(0) = 0
ami = pij(r).*log(pij(r)./pii(r)./pjj(r));
out = sum(ami);

end