% ------------------------------------------------------------------------------
% SY_LocalGlobal
% ------------------------------------------------------------------------------
% 
% Compares statistics measured in a local region of the time series to that
% measured of the full time series.
% 
%---INPUTS:
% y, the time series to analyze
% 
% lorp, the local subset of time series to study:
%             (i) 'l': the first n points in a time series,
%             (ii) 'p': an initial proportion of the full time series, n
%             (iii) 'unicg': n evenly-spaced points throughout the time series
%             (iv) 'randcg': n randomly-chosen points from the time series
%                               (chosen with replacement)
% 
% n, the parameter for the method specified above
% 
% 
%---OUTPUTS: the mean, standard deviation, median, interquartile range,
% skewness, kurtosis, AC(1), and SampEn(1,0.1).
% 
% This is not the most reliable or systematic operation because only a single
% sample is taken from the time series and compared to the full time series.
% A better approach would be to repeat over many local subsets and compare the
% statistics of these local regions to the full time series.
% 
% 
%---HISTORY:
% Ben Fulcher, September 2009
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = SY_LocalGlobal(y,lorp,n)

% ------------------------------------------------------------------------------
% Preliminaries
% ------------------------------------------------------------------------------

% Check input time series is z-scored:
if ~BF_iszscored(y)
    warning('The input time series should be z-scored (mean=%g and std-1=%g)', ...
                                            mean(y),std(y)-1)
end

% Set default l
if nargin < 2 || isempty(lorp)
    lorp = 'l';
end

% Set default n
if nargin < 3 || isempty(n)
    switch lorp
    case {'l','unicg','randcg'}
        n = 100; % 100 samples
    case 'p'
        n = 0.1; % 10% of the time series
    end
end

N = length(y); % Length of the time series

% ------------------------------------------------------------------------------
% Determine subset range to use: r
% ------------------------------------------------------------------------------
switch lorp
    case 'l'
        % Takes first n points of time series:
        r = (1:min(n,N));
    case 'p'
        % Takes initial proportion n of time series
    	r = (1:round(N*n));
    case 'unicg'
        % Takes n uniformly distributed points in time series:
        r = round(linspace(1,N,n));
    case 'randcg'
        % Reset the seed to the default value (for reproducibility):
        BF_ResetSeed('default');
        % Takes n random points in time series; there could be repeats:
        r = randi(N,n,1);
        % This is not very robust, as it's taking just a single stochastic
        % sample with a (possibly) large variance
    otherwise
        error('Unknown specifier, ''%s''',lorp);
end

% ------------------------------------------------------------------------------
% Compare statistics of this subset to those obtained from the full time series
% ------------------------------------------------------------------------------
out.absmean = abs(mean(y(r))); %/mean(y); % ** INPUT Y MUST BE Z-SCORED;
out.std = std(y(r)); %/std(y); % ** INPUT Y MUST BE Z-SCORED
out.median = median(y(r)); %/median(y); % if median is very small;; could be very noisy
out.iqr = abs(1-iqr(y(r)) / iqr(y));
out.skewness = abs(1-skewness(y(r)) / skewness(y)); % how far from true
out.kurtosis = abs(1-kurtosis(y(r)) / kurtosis(y)); % how far from true
out.ac1 = abs(1 - CO_AutoCorr(y(r),1) / CO_AutoCorr(y,1)); % how far from true
out.sampen101 = PN_sampenc(y(r),1,0.1,1) / PN_sampenc(y,1,0.1,1);


end