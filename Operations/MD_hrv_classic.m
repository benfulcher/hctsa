function out = MD_hrv_classic(y)
% Packages up classic HRV operations that calculate classical HRV 
% analysis measures from a given NN/RR time series in units of seconds.
% Adapted from code emailed from Max Little on 26/1/2009. 
% Implemented on 12/4/2010. 
% (Yep, a typically efficient turn-around time of over one year)
% Ben Fulcher, 12/4/2010.

% Standard defaults
diffy = diff(y);
N = length(y);

%% Calculate pNNx percentage
% pNNx: recommendation as per Mietus et. al. 2002, "The pNNx files: ...", Heart
% strange to do this for a z-scored time series...
% pnntime = 20;
Dy = abs(diffy) * 1000;
% exceed  = sum(Dy > pnntime);
out.pnn5  = sum(Dy > 5)/(N-1); % proportion of difference magnitudes greater than 0.005*sigma
out.pnn10 = sum(Dy > 10)/(N-1);
out.pnn20 = sum(Dy > 20)/(N-1);
out.pnn30 = sum(Dy > 30)/(N-1);
out.pnn40 = sum(Dy > 40)/(N-1);

% Calculate PSD
% [Pxx, F] = psd(series,1024,1,hanning(1024),512);
[Pxx, F] = periodogram(y,hann(N)); % periodogram with hanning window

% Calculate spectral measures such as subband spectral power percentage, LF/HF ratio etc.
% LF/HF: as per Malik et. al. 1996, "Heart Rate Variability"
LF_lo = 0.04; % /pi -- fraction of total power (max F is pi)
LF_hi = 0.15;
HF_lo = 0.15; 
HF_hi = 0.4;

fbinsize = F(2) - F(1);
indl  = ((F >= LF_lo) & (F <= LF_hi));
indh  = ((F >= HF_lo) & (F <= HF_hi));
indv  = (F <= LF_lo);
lfp   = fbinsize * sum(Pxx(indl));
hfp   = fbinsize * sum(Pxx(indh));
vlfp  = fbinsize * sum(Pxx(indv));
out.lfhf  = lfp / hfp;
total = fbinsize * sum(Pxx);
out.vlf   = vlfp/total * 100;
out.lf    = lfp/total * 100;
out.hf    = hfp/total * 100;

% Triangular histogram index
out.tri = length(y)/max(hist(y));

% Poincare plot measures: see
% "Do Existing Measures ... ", Brennan et. al. (2001), IEEE Trans Biomed Eng 48(11)
rmssd = std(diffy); % std of differenced series
sigma = std(y); % should be 1 for zscored time series
out.SD1 = 1/sqrt(2) * rmssd * 1000;
out.SD2 = sqrt(2 * sigma^2 - (1/2) * rmssd^2) * 1000;



end