% CO_BinaryStretch
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

function out = CO_BinaryStretch(x,stretchwhat)
% Ben Fulcher, 2009

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