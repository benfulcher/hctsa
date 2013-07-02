function out = TSTL_acf(y)
% finds autocorrelation function using fft of length l
% Uses code from the TSTOOL package
% Compares to the autocorrelation function computed directly
% using my code
% Ben Fulcher October 2009

% first 50 autocorrelations
try
    co_fft = data(acf(signal(y),100));
catch me
    if strcmp(me.message,'Matrix dimensions must agree.')
        out = NaN; return
    end
end

n = length(co_fft);
co_ben = zeros(n,1);
for i = 1:n
    co_ben(i) = CO_autocorr(y,i-1);
end
% keyboard

out = norm(co_ben-co_fft);

end