function p = MD_polvar(x, d, D)
% Calculate POLVARd measure for a time series of RR intervals x in units of seconds.
% d [USED TO BE(!)] the symbolic coding time difference in milliseconds.
% p = polvar(x, d)
% x - input time series
% d - symbolic coding (amplitude) difference
% p - probability of obtaining a sequence of consecutive ones/zeros
% D - word length (classically words of length 6)
% generally a binary coding

% I think first mention is in Wessel et al., PRE (2000), called Plvar
% Code obtained from Max Little on January 2009.
% Repackaged for more general applications to time series by Ben Fulcher, April 2010

% coding time difference, d
if nargin < 2 || isempty(d)
    d = 1;
end

% Word length, D
if nargin < 3 || isempty(D)
    D = 6;   % Look for symbol words of this length
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