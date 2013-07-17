function out = CO_tc3(y,tau)
% Implements the TC3 function
% This function is mentioned in the documentation for the TSTOOL package
% (full package available here: http://www.physik3.gwdg.de/tstool/)
% Input time series, y, and time lag, tau
% Ben Fulcher 15/11/2009

%% Set defaults:
if nargin < 2 || isempty(tau)
    tau = 'ac';
end

% Can set the time lag, tau, to be 'ac' or 'mi'
if strcmp(tau,'ac')
    tau = CO_fzcac(y);
    % tau is first zero crossing of the autocorrelation function
elseif strcmp(tau,'mi')
    tau = CO_firstmin(y,'mi');
    % tau is the first minimum of the automutual information function
end

yn = y(1:end-2*tau);
yn1 = y(1+tau:end-tau); % yn1, tau steps ahead
yn2 = y(1+2*tau:end); % yn2, 2*tau steps ahead

% The expression used in TSTOOL tc3:
out.raw = mean(yn.*yn1.*yn2)/abs(mean(yn.*yn1))^(3/2);

% The magnitude
out.abs = abs(out.raw);

% The numerator
out.num = mean(yn.*yn1.*yn2);
out.absnum = abs(out.num);

% The denominator
out.denom = abs(mean(yn.*yn1))^(3/2);

end