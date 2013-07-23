% SY_StatAv
% 
% The StatAv measure is a simple mean-stationarity metric that divides
% the time series into non-overlapping subsegments, calculates the mean in each
% of these segments and returns the standard deviation of this set of means.
% 
% "Heart rate control in normal and aborted-SIDS infants", S. M. Pincus et al.
% Am J. Physiol. Regul. Integr. Comp. Physiol. 264(3) R638 (1993)
% 
% INPUTS:
% 
% y, the input time series
% 
% whattype, the type of StatAv to perform:
%           (i) 'seg': divide the time series into n segments
%           (ii) 'len': divide the time series into segments of length n
% 
% n, either the number of subsegments ('seg') or their length ('len')
% 

function out = SY_StatAv(y,whattype,n)
% Ben Fulcher, 2009
% Might be nicer to use the 'buffer' function for this...?

if nargin < 2
    whattype = 'seg'; % divide into n segments by default
end

if nargin < 3
    n = 5; % use 5 segments
end

N = length(y); % time-series length

switch whattype
case 'seg'
    % divide time series into n segments
    M = zeros(n,1);
    p = floor(N/n);% lose the last N mod n data points
    
    for j = 1:n
        M(j) = mean(y(p*(j-1)+1:p*j));
    end

case 'len'
    if N > 2*n
        pn = floor(N/n);
        M = zeros(pn,1);
        for j = 1:pn
            M(j) = mean(y((j-1)*n+1:j*n));
        end
    else
        fprintf(1,'This time series (N = %u) is too short for StatAv(%s,%u)\n',N,whattype,n)
        out = NaN; return
    end
otherwise
    error('Error evaluating StatAv of type ''%s'', please select either ''seg'' or ''len''',whattype)
end

s = std(y); % should be 1 (for a z-scored time-series input)
sdav = std(M);
out = sdav/s;

end