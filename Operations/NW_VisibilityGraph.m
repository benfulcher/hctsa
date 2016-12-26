function out = NW_VisibilityGraph(y,meth,maxL)
% NW_VisibilityGraph    Visibility graph analysis of a time series.
%
% Constructs a visibility graph of the time series and returns various
% statistics on the properties of the resulting network.
%
% cf.: "From time series to complex networks: The visibility graph"
% Lacasa, Lucas and Luque, Bartolo and Ballesteros, Fernando and Luque, Jordi
% and Nuno, Juan Carlos P. Natl. Acad. Sci. USA. 105(13) 4972 (2008)
%
% "Horizontal visibility graphs: Exact results for random time series"
% Luque, B. and Lacasa, L. and Ballesteros, F. and Luque, J.
% Phys. Rev. E. 80(4) 046103 (2009)
%
%---INPUTS:
% y, the time series (a column vector)
%
% meth, the method for constructing:
% 			(i) 'norm': the normal visibility definition
% 			(ii) 'horiz': uses only horizonatal lines to link nodes/datums
%
% maxL, the maximum number of samples to consider. Due to memory constraints,
%               only the first maxL (6000 by default) points of time series are
%               analyzed. Longer time series are reduced to their first maxL
%               samples.
%
%---OUTPUTS:
% Statistics on the degree distribution, including the mode, mean,
% spread, histogram entropy, and fits to gaussian, exponential, and powerlaw
% distributions.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%% Preliminaries, check inputs
% ------------------------------------------------------------------------------
N = length(y); % time-series length

if size(y,2) > size(y,1), y = y'; end % make sure a column vector
if nargin < 2
    % compute the horizontal visibility graph by default
    meth = 'horiz';
end
if nargin < 3
    maxL = 5000; % crops time series longer than this maximum length
end

if N > maxL % too long to store in memory
    % ++BF changed on 8/3/2010 to reduce down to first maxL samples. In future,
    % could alter to take different subsets, or set a maximum distance range
    % allowed to make a link (using sparse), etc.
	warning(sprintf(['Time series (%u > %u) is too long for visibility graph...' ...
                ' Analyzing the first %u samples'],N,maxL,maxL));
    y = y(1:maxL);
    N = length(y); % new time-series length
end


y = y - min(y); % adjust so that minimum of y is at zero

% ------------------------------------------------------------------------------
%% Compute the visibility graph:
% ------------------------------------------------------------------------------
switch meth
	case 'norm'
        % Normal visibility graph:
        A = EZ_VisibilityGraph(y);

	case 'horiz'
        % Horizontal visibility graph

        A = zeros(N); % adjacency matrix
        yr = flipud(y); % reversed order time series

		for i = 1:N
			% Look forward to first blocker, then stop
            if i < N
    			nAhead = find(y(i+1:end) > y(i),1,'first');
                A(i,i+nAhead) = 1;
            end

            % Look back to the first hit, then stop
            if i > 1
    			nBack = find(yr(N-i+2:end) > yr(N-i+1),1,'first');
                A(i-nBack,i) = 1;
            end
		end

        % Symmetrize A:
        A = symmetrize(A);
    otherwise
        error('Unknown visibility graph method ''%s''',meth);
end

% ------------------------------------------------------------------------------
%%% Statistics on the output
% ------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%% Degree distribution: basic statistics
%-------------------------------------------------------------------------------
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

% ------------------------------------------------------------------------------
%% Fit distributions to degree distribution
% ------------------------------------------------------------------------------
% (1) Gauss1: Gaussian fit to degree distribution
try
    dgaussout = DN_SimpleFit(k,'gauss1',range(k)); % range(k)-bin single gaussian fit
catch emsg
    warning(sprintf('Error fitting gaussian distribution to data:\n%s',emsg.message))
    dgaussout = NaN;
end

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

% (2) Exponential1: Exponential fit to degree distribution
try
    dexpout = DN_SimpleFit(k,'exp1',range(k)); % range(k)-bin single exponential fit
catch emsg
    warning(sprintf('Error fitting exponential distribution to data:\n%s',emsg.message))
    dexpout = NaN;
end
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

% (3) Power1: Power-law fit to degree distribution
try
    dpowerout = DN_SimpleFit(k,'power1',range(k)); % range(k)-bin single power law fit
catch emsg
    warning(sprintf('Error fitting power-law distribution to data:\n%s',emsg.message))
    dpowerout = NaN;
end
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

% ------------------------------------------------------------------------------
%% Using likelihood now:
% ------------------------------------------------------------------------------
% Gaussian
out.gaussnlogL = normlike([mean(k),std(k)],k);

% Exp
out.expnlogL = explike(mean(k),k);

% Extreme Value Distribution
paramhat = evfit(k);
out.evparm1 = paramhat(1);
out.evparm2 = paramhat(2);
out.evnlogL = evlike(paramhat,k);

% ------------------------------------------------------------------------------
%% Entropy of distribution:
% ------------------------------------------------------------------------------
out.entropy = EN_DistributionEntropy(k,'hist','sqrt');

% Autocorrelations:
out.kac1 = CO_AutoCorr(k,1,'Fourier');
out.kac2 = CO_AutoCorr(k,2,'Fourier');
out.kac3 = CO_AutoCorr(k,3,'Fourier');
out.ktau = CO_FirstZero(k,'ac');

%-------------------------------------------------------------------------------
function A = symmetrize(A)
    % Symmetrize an upper triangular matrix:
    At = A';
    lowerT = logical(tril(ones(size(A))));
    A(lowerT) = At(lowerT);
end

function ind = findFirst(vector,threshold)
    % Find index of the first time a vector exceeds a threshold
    % -- not used because just as fast to use find(x,1,'first')
    for k = 1:length(vector)
        if vector(k)>threshold
            ind = k;
            return;
        end
    end
    ind = length(vector);
end

end
