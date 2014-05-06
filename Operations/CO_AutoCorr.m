% ------------------------------------------------------------------------------
% CO_AutoCorr
% ------------------------------------------------------------------------------
% 
% Computes the autocorrelation of an input time series, y, at a time-lag, tau
% 
%---INPUTS:
% y, a scalar time series column vector.
% tau, the time-delay. If tau is a scalar, returns autocorrelation for y at that
%       lag. If tau is a vector, returns autocorrelations for y at that set of
%       lags.
% WhatMethod, the method of computing the autocorrelation: 'Fourier',
%             'TimeDomainStat', or 'TimeDomain'.
%       
%---OUTPUT: the autocorrelation at the given time-lag.
%
%---HISTORY:
% Ben Fulcher, 2014-03-24. Added multiple definitions for computing the
%       autocorrelation Computing mean/std across the full time series makes a
%       significant difference for short time series, but can produce values
%       outside [-1,+1]. The filtering-based method used by Matlab's autocorr,
%       is probably the best for short time series, and is now implemented here.
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

function out = CO_AutoCorr(y,tau,WhatMethod)

% ------------------------------------------------------------------------------
% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(tau)
    tau = 1; % Use a lag of 1 by default
end

if nargin < 3 || isempty(WhatMethod)
    WhatMethod = 'TimeDomain';
end

% ------------------------------------------------------------------------------
% Evaluate the time-series autocorrelation
% ------------------------------------------------------------------------------

switch WhatMethod
case 'Fourier'
    % ------------------------------------------------------------------------------
    % Estimation based on Matlab function autocorr, based on method of:
    % [1] Box, G. E. P., G. M. Jenkins, and G. C. Reinsel. Time Series
    %     Analysis: Forecasting and Control. 3rd edition. Upper Saddle River,
    %     NJ: Prentice-Hall, 1994.

    nFFT = 2^(nextpow2(length(y))+1);
    F = fft(y-mean(y),nFFT);
    F = F.*conj(F);
    acf = ifft(F);
    acf = acf./acf(1); % Normalize
    acf = real(acf);
    
    % acf = acf(2:end);

    if max(tau) > length(acf)-1 % -1 because acf(1) is lag 0
        warning(1,'Time lag %u is too long for time-series length %u',max(tau),length(y))
    end
    if min(tau) < 0
        warning(1,'Negative time lags not applicable')
    end
    
    out = zeros(length(tau),1);
    for i = 1:length(tau)
        if (tau(i) > length(acf)-1) || (tau(i) < 0)
            out(i) = NaN;
        else
            out(i) = acf(tau(i)+1);
        end
    end


case 'TimeDomainStat'
    % ------------------------------------------------------------------------------
    % Assume a stationary process and estimate mean and standard deviation
    % from the full time series, but then values can be outside [-1,1]...
    
    N = length(y); % time-series length
    sigma2 = var(y); % time-series standard deviation
    mu = mean(y); % time-series mean

    % Define the time-domain autocorrelation function:
    ACFy = @(tau) mean((y(1:N-tau) - mu).*(y(tau+1:N) - mu))/sigma2;

    if length(tau) == 1
        % Output a single value at the given time-lag
        out = ACFy(tau);
    else
        % Output values over a range of time lags
        out = arrayfun(@(x)ACFy(x),tau);
    end
    
case 'TimeDomain'
    % ------------------------------------------------------------------------------
    % Estimate mean and standard deviation from each portion of the time series

    N = length(y); % time-series length

    if length(tau) == 1 % Output a single value at the given time-lag
        out = mean((y(1:N-tau) - mean(y(1:N-tau))).* ...
                (y(tau+1:N) - mean(y(tau+1:N))))/std(y(1:N-tau))/std(y(tau+1:N));
    else
        % Output values over a range of time lags
        out = zeros(length(tau),1);
        for i = 1:length(tau)
            out(i) = mean((y(1:N-tau(i)) - mean(y(1:N-tau(i)))).*(y(tau(i)+1:N) ...
                        - mean(y(tau(i)+1:N))))/std(y(1:N-tau(i)))/std(y(tau(i)+1:N));
        end
    end

otherwise
    error('Unknown autocorrelation estimation method ''%s''',WhatMethod)

end