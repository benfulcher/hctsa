% MD_hrv_classic
% 
% Packages up a bunch of classic heart rate variability (HRV) statistics and
% applies them to the input time series.
% 
% Assumes an NN/RR time series in units of seconds.
% 
% INPUTS:
% y, the input time series.
% 
% Includes:
%  (i) pNNx
%  cf. "The pNNx files: re-examining a widely used heart rate variability
%           measure", J.E. Mietus et al., Heart 88(4) 378 (2002)
% 
%  (ii) Power spectral density ratios in different frequency ranges
%   cf. "Heart rate variability: Standards of measurement, physiological
%       interpretation, and clinical use",
%       M. Malik et al., Eur. Heart J. 17(3) 354 (1996)
% 
%  (iii) Triangular histogram index, and
%  
%  (iv) Poincare plot measures
%  cf. "Do existing measures of Poincare plot geometry reflect nonlinear
%       features of heart rate variability?"
%       M. Brennan, et al., IEEE T. Bio.-Med. Eng. 48(11) 1342 (2001)
%  
% Code is heavily derived from that provided by Max A. Little:
% http://www.maxlittle.net/
% 
function out = MD_hrv_classic(y)
% Ben Fulcher, 12/4/2010.

% Standard defaults
diffy = diff(y);
N = length(y); % time-series length

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