function [mu, sigma, lfhf, pnnx, vlf, lf, hf, rmssd, tri, SD1, SD2] = hrv_classic(series)
% Calculates classical HRV analysis measures from a given NN/RR time series
% in units of seconds.
% [mu, sigma, lfhf, pnnx, vlf, lf, hf, rmssd, tri, SD1, SD2] = hrv_classic(series)

% Standard defaults

% LF/HF: as per Malik et. al. 1996, "Heart Rate Variability"
LF_lo = 0.04;
LF_hi = 0.15;
HF_lo = 0.15; 
HF_hi = 0.4;

% pNNx: recommendation as per Mietus et. al. 2002, "The pNNx files: ...", Heart
pnntime = 20;

% Calculate simple statistics
mu    = mean(series);
sigma = std(series);

diffseries = diff(series);

% Simple statistics of difference series
rmssd   = std(diffseries);

% Calculate pNNx percentage
dseries = abs(diffseries) * 1000;
exceed  = sum(dseries > pnntime);
pnnx    = exceed / (length(series)-1);

% Calculate PSD
[Pxx, F] = psd(series,1024,1,hanning(1024),512);

% Calculate spectral measures such as subband spectral power percentage, LF/HF ratio etc.
fbinsize = F(2) - F(1);
indl  = find((F >= LF_lo) & (F <= LF_hi));
indh  = find((F >= HF_lo) & (F <= HF_hi));
indv  = find(F <= LF_lo);
lfp   = fbinsize * sum(Pxx(indl));
hfp   = fbinsize * sum(Pxx(indh));
vlfp  = fbinsize * sum(Pxx(indv));
lfhf  = lfp / hfp;
total = fbinsize * sum(Pxx);
vlf   = vlfp/total * 100;
lf    = lfp/total * 100;
hf    = hfp/total * 100;

% Triangular histogram index
tri   = length(series)/max(hist(series));

% Poincare plot measures: see
% "Do Existing Measures ... ", Brennan et. al. (2001), IEEE Trans Biomed Eng 48(11)
SD1 = 1/sqrt(2) * rmssd * 1000;
SD2 = sqrt(2 * sigma^2 - (1/2) * rmssd^2) * 1000;
