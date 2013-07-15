function out = TSTL_delaytime(y,maxdelay,past)
% Uses TSTOOL code delaytime
% y: column vector of time series data
% maxdelay: maximum value of the delay to consider (can also specify a
%           proportion of time series length)
% past: the TSTOOL documentation describes this parameter as "?", which is
%       relatively uninformative.
% Uses the method of "Parlitz and Wichard", which is unreferenced. So I
% really have no idea what's going on with this algorithm!!
% I know it's a stochastic algorithm, so it must rely on some random
% sampling of the input time series... I'll just return some statistics
% Ben Fulcher, November 2009

%% Preliminaries
N = length(y); % length of time series
try
    s = signal(y); % convert to a signal for TSTOOL
catch
    error('Error converting time series to signal class using TSTOOL function ''signal''')
end

% (1) Maximum delay, maxdelay
if nargin < 2 || isempty(maxdelay)
    maxdelay = 0.2; % 1/5 the length of the time series
end
if maxdelay < 1 && maxdelay > 0
    maxdelay = round(N*maxdelay); % specify a proportion of time series length
end

if maxdelay < 10,
    maxdelay = 10;
    fprintf(1,'Maxdelay set to its minimum: delaytime = 10\n')
    % necessary for output statistics
end

%% Run
tau = data(delaytime(s,maxdelay,past));
% plot(tau);

%% Output Statistics
% tau tends to start low and then rise to some (noisy) value
out.tau1 = tau(1);
out.tau2 = tau(2);
out.tau3 = tau(3);
out.difftau12 = tau(2)-tau(1);
out.difftau13 = tau(3)-tau(1);
out.meantau = mean(tau);
out.stdtau = std(tau);
out.mintau = min(tau);
out.maxtau = max(tau);


end