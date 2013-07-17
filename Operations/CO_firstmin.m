function out = CO_firstmin(y,minwhat)
% (previously CO_fmmi)
% (1) Finds the first minimum of the automutual information function
% Uses Rudy Moddemeijer's information.m code (RM_information.m here) that may or may not be great...
% (2) Outputs the first minimum of the linear autocorrelation function using CO_autocorr
% Clumsily coded by Ben Fulcher, 2008

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
end

% Go through lags until find a minimum
a = zeros(N-1,1); % preallocate space for the autocorrelation function
for i = 0:N-1
    a(i+1) = corrfn(i); % calculate the auto-correlation at this lag
            
    if (i > 1) && (a(i-1)-a(i) > 0) && (a(i)-a(i+1) < 0); % minimum
        out = i-1; % I found the first minimum!
        return
    end
end

% If there is no minimum in the automutual information
out = N;

end

