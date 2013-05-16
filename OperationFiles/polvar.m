function p = polvar(x, d, D)
% Calculate POLVARd measure for a time series of RR intervals x in units of seconds.
% d is the symbolic coding time difference in milliseconds.
% p = polvar(x, d)
% x - input time series
% d - symbolic coding time difference in milliseconds
% p - probability of obtaining a sequence of consecutive ones/zeros
dx   = abs(diff(x));
xsym = (dx >= (d / 1000));
N    = length(xsym);
% D    = 6;
% Look for symbol words of length D
zseq = zeros(D, 1);
oseq = ones(D, 1);

% Search for D consecutive zeros/ones
i      = 1;
pc     = 0;
%seqcnt = 0;
while (i <= (N - D))
    xseq = xsym(i:(i+D-1));
%    seqcnt = seqcnt + 1;
    if ((sum(xseq == zseq) == D) | (sum(xseq == oseq) == D))
        pc = pc + 1;
        i  = i + D;
    else
        i = i + 1;
    end
end

p = pc / length(xsym);
%p = pc / seqcnt;
