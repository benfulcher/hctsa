% CO_FirstZero
% 
% Returns the first zero-crossing of a given autocorrelation function.
% 
% y, the input time series
% corrfn, the self-correlation function to measure:
%         (i) 'ac': normal linear autocorrelation function. Uses CO_AutoCorr to
%                   calculate autocorrelations.
% maxtau, a maximum time-delay to search up to
% 
% In future, could add an option to return the point at which the function
% crosses the axis, rather than the first integer lag at which it has already
% crossed (what is currently implemented)
% 

function out = CO_FirstZero(y,corrfn,maxtau)
% Ben Fulcher, 2008

N = length(y); % the length of the time series

if nargin < 2 || isempty(corrfn)
    corrfn = 'ac'; % autocorrelation by default
end
if nargin < 3 || isempty(maxtau)
    maxtau = N; % search up to a maximum of the length of the time series
    % maxtau = 400; % searches up to this maximum time lag
    % maxtau = min(maxtau,N); % searched up to the length of the time series if this is less than maxtau
end

% Select the self-correlation function as an inline function
% Eventually could add additional self-correlation functions
switch corrfn
case 'ac'
    Corr_fn = @(x) CO_AutoCorr(y,x); % autocorrelation at time lag x
otherwise
    error('Unknown correlation function ''%s''',corrfn);
end

% Calculate autocorrelation at increasing lags, until you find a negative one
for tau = 1:maxtau-1
    if Corr_fn(tau) < 0 % we know it starts positive (1), so first negative will be the zero-crossing
        out = tau; return
    end
end

% If haven't left yet, set output to maxtau
out = maxtau;

end