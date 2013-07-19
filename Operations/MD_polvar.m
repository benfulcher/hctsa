% MD_polvar
% 
% Calculates the POLVARd measure for a time series.
% The first mention may be in Wessel et al., PRE (2000), called Plvar
% cf. "Short-term forecasting of life-threatening cardiac arrhythmias based on
% symbolic dynamics and finite-time growth rates",
%       N. Wessel et al., Phys. Rev. E 61(1) 733 (2000)
% 
% The output from this function is the probability of obtaining a sequence of
% consecutive ones or zeros.
% 
% Although the original measure used raw thresholds, d, on RR interval sequences
% (measured in milliseconds), this code can be applied to general z-scored time
% series. So now d is not the time difference in milliseconds, but in units of
% std.
% 
% The measure was originally applied to sequences of RR intervals and this code
% is heavily derived from that provided by Max A. Little, January 2009.
% cf. http://www.maxlittle.net/
% 
% INPUTS:
% x, the input time series
% d, the symbolic coding (amplitude) difference,
% D, the word length (classically words of length 6).
% 
% OUPUT:
% p - probability of obtaining a sequence of consecutive ones/zeros

function p = MD_polvar(x, d, D)
% Ben Fulcher, April 2010

% coding time difference, d
if nargin < 2 || isempty(d)
    d = 1;
end

% Word length, D
if nargin < 3 || isempty(D)
    D = 6;   % Look for symbol words of this default length
end

dx = abs(diff(x)); % absolute differences in consecutive values of the time series
N = length(dx); % number of differences in the input time series

% binary representation of time series based on consecutive changes being
% greater than d/1000...
xsym = (dx >= d); % consecutive differences exceed than some threshold, d
zseq = zeros(D,1);
oseq = ones(D,1);

% Search for D consecutive zeros/ones
i  = 1;
pc = 0;

%seqcnt = 0;
while (i <= (N - D))
    xseq = xsym(i:(i+D-1));
%    seqcnt = seqcnt + 1;
    if (sum(xseq == zseq) == D) || (sum(xseq == oseq) == D)
        pc = pc + 1;
        i  = i + D;
    else
        i = i + 1;
    end
end

p = pc / N;
%p = pc / seqcnt;

end