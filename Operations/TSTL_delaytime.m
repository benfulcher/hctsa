% TSTL_delaytime
% 
% Uses the TSTOOL code delaytime, that computes an optimal delay time using the
% method of Parlitz and Wichard (this method is specified in the TSTOOL
% documentation but without reference).
% 
% TSTOOL: http://www.physik3.gwdg.de/tstool/
% 
% INPUTS:
% y, column vector of time series data
% 
% maxdelay, maximum value of the delay to consider (can also specify a
%           proportion of time series length)
%           
% past, the TSTOOL documentation describes this parameter as "?", which is
%       relatively uninformative.
% 
% 
% It's a stochastic algorithm, so it must rely on some random sampling of the
% input time series... A bit of a strange one, but I'll return some statistics
% and see what they do.
% 

function out = TSTL_delaytime(y,maxdelay,past)
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