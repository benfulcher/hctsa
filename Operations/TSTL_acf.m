function out = TSTL_acf(y)
% finds autocorrelation function using fft of length l
% Uses code from the TSTOOL package
% Compares to the autocorrelation function computed directly
% using my code
% Ben Fulcher October 2009

% first 50 autocorrelations
% try
co_fft = data(acf(signal(y),100));
% catch emsg
%     if strcmp(emsg.message,'Matrix dimensions must agree.')
%         out = NaN; return
%     end
% end

nlags = length(co_fft);
co_ben = zeros(nlags,1);
for i = 1:nlags
    co_ben(i) = CO_autocorr(y,i-1);
end

out = norm(co_ben - co_fft);

end