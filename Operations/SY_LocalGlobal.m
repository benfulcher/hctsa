% SY_LocalGlobal
% 
% Compares statistics measured in a local region of the time series to that
% measured of the full time series.
% 
% INPUTS:
% y, the time series to analyze
% 
% lorp, the local subset of time series to study:
%             (i) 'l': the first n points in a time series,
%             (ii) 'p': an initial proportion of the full time series, n
%             (iii) 'unicg': n evenly-spaced points throughout the time series
%             (iv) 'randcg': n randomly-chosen points from the time series (chosen with replacement)
% 
% n, the parameter for the method specified above
% 
% Statistics are outputted in a structure: the mean, standard deviation, median,
% interquartile range, skewness, kurtosis, AC(1), and SampEn(1,0.1).
% 
% This is not the most reliable or systematic operation because only a single
% sample is taken from the time series and compared to the full time series.
% A better approach would be to repeat over many local subsets and compare the
% statistics of these local regions to the full time series.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
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
% Ben Fulcher, September 2009

% Check z-scored time series
if ~BF_iszscored(y)
    warning('The input time series should be z-scored')
end

if nargin < 2 || isempty(lorp)
    lorp = 'l';
end

if nargin < 3 || isempty(n)
    switch lorp
    case {'l','unicg','randcg'}
        n = 100; % 100 samples
    case 'p'
        n = 0.1; % 10% of the time series
    end
end
N = length(y); % length of the time series

% Determine subset range to use: r
switch lorp
    case 'l'
        r = (1:min(n,N)); % takes first n points of time series
    case 'p'
    	r = (1:round(N*n)); % takes initial proportion n of time series
    case 'unicg'
        r = round(linspace(1,N,n)); % takes n uniformly distributed points in time series
    case 'randcg'
        r = randi(N,n,1); % takes n random points in time series; there could be repeats
        % This is quite unrobust, as it's taking just a single sample from
        % a test with a (possibly) large variance
    otherwise
        error('Unknown specifier, ''%s''',lorp);
end


% Compare this subset to the full value
out.mean = abs(mean(y(r))); %/mean(y); % Y SHOULD BE Z-SCORED;;
out.std = std(y(r)); %/std(y); % Y SHOULD BE Z-SCORED;;
out.median = median(y(r)); %/median(y); % if median is very small;; could be very noisy
out.iqr = abs(1-iqr(y(r)) / iqr(y));
out.skewness = abs(1-skewness(y(r)) / skewness(y)); % how far from true
out.kurtosis = abs(1-kurtosis(y(r)) / kurtosis(y)); % how far from true
out.ac1 = abs(1-CO_AutoCorr(y(r),1) / CO_AutoCorr(y,1)); % how far from true
out.sampen101 = PN_sampenc(y(r),1,0.1,1) / PN_sampenc(y,1,0.1,1);


% switch wing
%     case 'mean'
%         out = abs(mean(y(r))); %/mean(y); % Y SHOULD BE Z-SCORED;;
%     case 'std'
%         out = std(y(r)); %/std(y); % Y SHOULD BE Z-SCORED;;
%     case 'median'
%         out = median(y(r)); %/median(y); % if median is very small;; could be very noisy
%     case 'iqr'
%         out = abs(1-iqr(y(r)) / iqr(y));
%     case 'skewness'
%         out = abs(1-skewness(y(r)) / skewness(y)); % how far from true
%     case 'kurtosis'
%         out = abs(1-kurtosis(y(r)) / kurtosis(y)); % how far from true
%     case 'AC1'
%         out = abs(1-CO_AutoCorr(y(r),1) / CO_AutoCorr(y,1)); % how far from true
%     case 'SampEn1_01' % computationally expensive to calculate this full one each time...
%         out = PN_sampenc(y(r),1,0.1,1) / PN_sampenc(y,1,0.1,1);
%     otherwise
%         error('Unknwon statistic ''%s''',wing);
% end



end