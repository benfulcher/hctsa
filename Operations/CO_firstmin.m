% CO_FirstMin
% 
% Returns the time at which the first minimum in a given correlation function
% occurs.
% 
% INPUTS:
% y, the input time series
% minwhat, the type of correlation to minimize: either 'ac' for autocorrelation,
%           or 'mi' for automutual information
% 
% Note that selecting 'ac' is unusual operation: standard operations are the
% first zero-crossing of the autocorrelation (as in CO_fzcac), or the first
% minimum of the mutual information function ('mi').
%
% Uses Rudy Moddemeijer's RM_information.m code that may or may not be great...
% 

function out = CO_FirstMin(y,minwhat)
% Ben Fulcher, 2008


if nargin < 2 || isempty(minwhat)
    minwhat = 'mi'; % mutual information
end

N = length(y); % time-series length

switch minwhat
case 'mi'
    % automutual information implemented as RM_information
    corrfn = @(x) RM_information(y(1:end-x), y(1+x:end));
case 'ac'
    % autocorrelation implemented as CO_autocorr
    corrfn = @(x) CO_autocorr(y,x);
otherwise
    error('Unknown correlation function ''%s''',minwhat);
end

% Go through lags until find a minimum
a = zeros(N-1,1); % preallocate space for the autocorrelation function
for i = 0:N-1
    a(i+1) = corrfn(i); % calculate the auto-correlation at this lag
            
    if (i > 1) && (a(i-1) - a(i) > 0) && (a(i) - a(i+1) < 0); % minimum
        out = i-1; % I found the first minimum!
        return
    end
end

% If there is no minimum in the automutual information
out = N;

end

