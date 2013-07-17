function out = CO_trev(y,tau)
% Implements the trev function
% (which is mentioned in the TSTOOL package, and elsewhere...)
% (the TSTOOL package is available here: http://www.physik3.gwdg.de/tstool/)
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

yn = y(1:end-tau);
yn1 = y(1+tau:end); % yn, tau steps ahead

% The expression used in TSTOOL
out.raw = mean((yn1-yn).^3)/(mean((yn1-yn).^2))^(3/2);

% The magnitude
out.abs = abs(out.raw);

% The numerator
out.num = mean((yn1-yn).^3);
out.absnum = abs(out.num);

% The denominator
out.denom = (mean((yn1-yn).^2))^(3/2);

end