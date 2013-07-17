function out = TSTL_acf(y)
% Calculates the autocorrelation function using fft of length l
% using code from the TSTOOL package then compares the results
% to the autocorrelation function computed directly using my
% own code.
% Ben Fulcher October 2009

% first 50 autocorrelations
co_fft = data(acf(signal(y),100));

nlags = length(co_fft);
co_ben = zeros(nlags,1);
for i = 1:nlags
    co_ben(i) = CO_autocorr(y,i-1);
end

out = norm(co_ben - co_fft);

end