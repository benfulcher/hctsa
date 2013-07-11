function out = NW_visgraphben(y,meth)
% Constructs a visibility graph from the time series (c.f., Lacasa 2008, PNAS)
% Returns various statistics on the graph.
% INPUTS:
% y: the time series (a column vector)
% meth: the method for constructing:
% 			(a) 'norm': the normal visibility definition
% 			(b) 'horiz': uses only horizonatal lines to link nodes/datums (Luque et al. PRE 2009)
% Ben Fulcher October 2009
% Ben Fulcher 5/6/2010 -- changed some sparse/full settings

%% Preliminaries
N = length(y); % length of time series

maxL = 6000;
if N > maxL % too long
    % ++BF changed on 8/3/2010 to reduce down to first maxL samples. In future,
    % could alter to take different subsets, or set a maximum distance range
    % allowed to make a link (using sparse), etc.
	fprintf(1,'Time series (%u > %u) is too long for visibility graph... Analyzing the first %u samples',N,maxL,maxL);
    y = y(1:maxL);
    N = length(y); % new time-series length
% 	return	
end
% if N > 4000 % needs too much memory -- need to use (slower to index) sparse representation
% %     A = sparse(N,N); % sparse is very slow for matricies that end up having
% %     many zeros
%     A = zeros(N); % adjacency matrix -- faster for shorter time series
% else
%     A = zeros(N); % adjacency matrix -- faster for shorter time series
% end
% % end

A = zeros(N); % adjacency matrix
y = y - min(y); % adjust so that minimum of y is at zero (no real reason other than aesthetics ;-))

%% Calculate the visibility graph
switch meth
	case 'norm'
		for i = 1 : N-1
			% compute all subsequent gradients
			deltay = y(i+1:end) - ones(N-i,1)*y(i); % vector of deltay's
			deltat = (1:N-i)'; % time from current reference i
			m = deltay./deltat; % gradients
			cummax = zeros(N-i,1);
			for j = 1:N-i
				cummax(j) = max(m(1:j));
			end
			links = (m >= cummax);
	
			% so we have 'links' which are the time series points following it in the series
			% that are visible from it (1 means the following point, and should always be included)
			% Store this information in the adjacency matrix, A
			A(i,i+1:end) = links';
		end
	case 'horiz'
		for i = 1:N-1
			% find first blocker and stop
			stophere = find(y(i+1:end) > y(i),1,'first');
			links = zeros(N-i,1);
			links(1:stophere) = 1;
			% so we have 'links' which are the time series points following it in the series
			% that are visible from it (1 means the following point, and should always be included)
			% Store this information in the adjacency matrix, A
			A(i,i+1:end) = links';
		end
end

% symmetrize A crudely: (works since lower triangle is zeros)
if N <= 5000, A = sparse(A); end
A = A + A'; A(A > 0) = 1;

%%% Statistics on the output
% spy(A);
%% Degree distribution: basic statistics
k = sum(A); % the degree distribution
k = full(k);

out.modek = mode(k); % mode of degree distribution
out.propmode = sum(k == mode(k))/sum(k);
out.meank = mean(k); % mean number of links per node
out.mediank = median(k); % median number of links per node
out.stdk = std(k); % std of k
out.maxk = max(k); % maximum degree
out.mink = min(k); % minimum degree
out.rangek = range(k); % range of degree distribution
out.iqrk = iqr(k); % interquartile range of degree distribution
out.maxonmedian = max(k)/median(k); % max on median (indicator of outlier)
out.ol90 = mean(k(k>=quantile(k,0.05) & k<=quantile(k,0.95)))/mean(k);
out.olu90 = (mean(k(k>=quantile(k,0.95)))-mean(k))/std(k); % top 5% of points are
                                                       % how far from mean (in std units)?

%% Fit distributions to degree distribution
% % (1) Gauss1: Gaussian fit to degree distribution
dgaussout = MF_M_mtlbfit(k,'gauss1',range(k)); % range(k)-bin single gaussian fit
if ~isstruct(dgaussout) && isnan(dgaussout)
    out.dgaussk_r2 = NaN;
    out.dgaussk_adjr2 = NaN;
    out.dgaussk_rmse = NaN;
    out.dgaussk_resAC1 = NaN;
    out.dgaussk_resAC2 = NaN;
    out.dgaussk_resruns = NaN;
else
    out.dgaussk_r2 = dgaussout.r2; % rsquared
    out.dgaussk_adjr2 = dgaussout.adjr2; % degrees of freedom-adjusted rsqured
    out.dgaussk_rmse = dgaussout.rmse;  % root mean square error
    out.dgaussk_resAC1 = dgaussout.resAC1; % autocorrelation of residuals at lag 1
    out.dgaussk_resAC2 = dgaussout.resAC2; % autocorrelation of residuals at lag 2
    out.dgaussk_resruns = dgaussout.resruns; % runs test on residuals -- outputs p-value
end

% % (2) Exponential1: Exponential fit to degree distribution
dexpout = MF_M_mtlbfit(k,'exp1',range(k)); % range(k)-bin single exponential fit
if ~isstruct(dexpout) && isnan(dexpout)
    out.dexpk_r2 = NaN;
    out.dexpk_adjr2 = NaN;
    out.dexpk_rmse = NaN;
    out.dexpk_resAC1 = NaN;
    out.dexpk_resAC2 = NaN;
    out.dexpk_resruns = NaN;
else
    out.dexpk_r2 = dexpout.r2; % rsquared
    out.dexpk_adjr2 = dexpout.adjr2; % degrees of freedom-adjusted rsqured
    out.dexpk_rmse = dexpout.rmse;  % root mean square error
    out.dexpk_resAC1 = dexpout.resAC1; % autocorrelation of residuals at lag 1
    out.dexpk_resAC2 = dexpout.resAC2; % autocorrelation of residuals at lag 2
    out.dexpk_resruns = dexpout.resruns; % runs test on residuals -- outputs p-value
end

% % (3) Power1: Power-law fit to degree distribution
dpowerout = MF_M_mtlbfit(k,'power1',range(k)); % range(k)-bin single power law fit
if ~isstruct(dpowerout) && isnan(dpowerout)
	out.dpowerk_r2 = NaN;
	out.dpowerk_adjr2 = NaN;
	out.dpowerk_rmse = NaN;
	out.dpowerk_resAC1 = NaN;
	out.dpowerk_resAC2 = NaN;
	out.dpowerk_resruns = NaN;
else
	out.dpowerk_r2 = dpowerout.r2; % rsquared
	out.dpowerk_adjr2 = dpowerout.adjr2; % degrees of freedom-adjusted rsqured
	out.dpowerk_rmse = dpowerout.rmse;  % root mean square error
	out.dpowerk_resAC1 = dpowerout.resAC1; % autocorrelation of residuals at lag 1
	out.dpowerk_resAC2 = dpowerout.resAC2; % autocorrelation of residuals at lag 2
	out.dpowerk_resruns = dpowerout.resruns; % runs test on residuals -- outputs p-value
end


% Using likelihood now:
% (1) Gauss
[muhat, sigmahat] = normfit(k);
out.gaussmu = muhat;
out.gausssigma = sigmahat;
out.gaussnlogL = normlike([muhat,sigmahat],k);

% (2) Exp
muhat = expfit(k);
out.expmu = muhat;
out.expnlogL = explike(muhat,k);

% (3) Poisson
lambdahat = poissfit(k);
out.explambda = lambdahat;

% (4) Extreme Value Distribution
paramhat = evfit(k);
out.evparm1 = paramhat(1);
out.evparm2 = paramhat(2);
out.evnlogL = evlike(paramhat,k);

%% Entropy of distribution
% for a range of bins

binr = (10:100); % range of nbins to try
h = zeros(size(binr));
for i = 1:length(binr);
    [n, x] = hist(k,binr(i));
    n = n/N;
    h(i) = - sum(n(n>0).*log(n(n>0)));
end
out.maxent = max(h);
out.minnbinmaxent = binr(find(h == max(h),1,'first'));
out.meanent = mean(h);
diffh = diff(h);
out.meanchent = mean(diffh(diffh~=0));

out.kac1 = CO_autocorr(k,1);
out.kac2 = CO_autocorr(k,2);
out.kac3 = CO_autocorr(k,3);
out.ktau = CO_fzcac(k);

end