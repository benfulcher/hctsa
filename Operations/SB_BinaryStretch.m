% SB_BinaryStretch
% 
% Measures the longest stretch of consecutive zeros or ones in a symbolized time
% series as a proportion of the time-series length.
% 
% The time series is symbolized to a binary string by whether it's above (1) or
% below (0) its mean.
% 
% It doesn't actually measure this correctly, due to an error in the code, but
% it's still kind of an interesting operation...?!
% 
% INPUTS:
% x, the input time series
% stretchwhat, (i) 'lseq1', measures something related to consecutive 1s
%              (ii) 'lseq0', measures something related to consecutive 0s
% 

function out = SB_BinaryStretch(x,stretchwhat)
% Ben Fulcher, 2009

if nargin < 2 || isempty(stretchwhat)
    stretchwhat = 'lseq1'; % by default
end

N = length(x); % length of the time series
x(x > 0) = 1;
x(x <= 0) = 0;

switch stretchwhat
    case 'lseq1'
        % longest stretch of 1s (this doesn't actually measure this!)
        out = max(diff(BF_sgnchange(diff(find(x == 1))-1.5,1)))/N;
    case 'lseq0'
        % longest stretch of 0s (this doesn't actualy measure this!)
        out = max(diff(BF_sgnchange(diff(find(x == 0))-1.5,1)))/N;
    otherwise
        error('Unknown input %s',stretchwhat)
end

if isempty(out)
    out = 0;
end

end