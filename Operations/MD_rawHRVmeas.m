function out = MD_rawHRVmeas(x)
% Ben Fulcher 24/2/2011 -- want to do this on raw RR intervals; may not
% make sense to other time series
% Note that pNNx is done in MD_pNN.m
% This code is very much derived from Max Little's hrv_classic.m code

N = length(x);

% Triangular histogram index
out.tri10 = N/max(hist(x,10));
out.tri20 = N/max(hist(x,20));
out.trisqrt = N/max(hist(x,sqrt(N)));


% Poincare plot measures: see
% "Do Existing Measures ... ", Brennan et. al. (2001), IEEE Trans Biomed Eng 48(11)
diffx = diff(x);
out.SD1 = 1/sqrt(2) * std(diffx) * 1000;
out.SD2 = sqrt(2 * var(x) - (1/2) * std(diffx)^2) * 1000;


% % Calculate PSD
% Hs = spectrum.periodogram('Hann');
% Hpsd = psd(Hs,x);
% keyboard
% % [Pxx, F] = spectrum.psd(series,1024,1,hanning(1024),512);
% 
% % LF/HF: as per Malik et. al. 1996, "Heart Rate Variability"
% LF_lo = 0.04;
% LF_hi = 0.15;
% HF_lo = 0.15; 
% HF_hi = 0.4;
% 
% % Calculate spectral measures such as subband spectral power percentage, LF/HF ratio etc.
% fbinsize = F(2) - F(1);
% indl  = ((F >= LF_lo) & (F <= LF_hi));
% indh  = ((F >= HF_lo) & (F <= HF_hi));
% indv  = (F <= LF_lo);
% lfp   = fbinsize * sum(Pxx(indl));
% hfp   = fbinsize * sum(Pxx(indh));
% vlfp  = fbinsize * sum(Pxx(indv));
% lfhf  = lfp / hfp;
% total = fbinsize * sum(Pxx);
% vlf   = vlfp/total * 100;
% lf    = lfp/total * 100;
% hf    = hfp/total * 100;

end
