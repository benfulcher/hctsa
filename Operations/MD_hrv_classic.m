function out = MD_hrv_classic(y)
% MD_hrv_classic   Classic heart rate variability (HRV) statistics.
%
% Typically assumes an NN/RR time series in units of seconds. Returns the pNNx
% measures (the proportion of successive differences larger than x/1000), the
% proportions of power in the very-low, low and high frequency bands of a
% Hann-windowed periodogram and the ratio of low to high, the triangular
% histogram index, and Poincare plot measures. The frequency bands are applied
% to the periodogram's normalized frequency in radians per sample, not in hertz.
%
% ---INPUTS:
% y, the input time series
%
% ---OUTPUTS:
% pnn5, pnn10, pnn20, pnn30, pnn40: proportion of successive differences of y
%       larger than 0.005, 0.010, 0.020, 0.030 and 0.040 (x/1000)
% lfhf, the ratio of power in the low-frequency band (0.04 to 0.15) to that in
%       the high-frequency band (0.15 to 0.4)
% vlf, the percentage of total power in the very-low-frequency band (below
%       0.04)
% lf, the percentage of total power in the low-frequency band (0.04 to 0.15)
% hf, the percentage of total power in the high-frequency band (0.15 to 0.4)
% tri, the triangular index: the length of the series divided by the count in
%       the fullest of 10 equal-width histogram bins
% SD1, the short-term variability from the Poincare plot: the standard
%       deviation of successive differences, divided by sqrt(2), times 1000
% SD2, the long-term variability from the Poincare plot,
%       sqrt(2*std(y)^2 - std(diff(y))^2/2), times 1000
%
% ---REFERENCES:
% Mietus et al., "The pNNx files: re-examining a widely used heart rate
% variability measure", Heart 88(4), 378 (2002).
%
% Malik et al., "Heart rate variability: Standards of measurement,
% physiological interpretation, and clinical use", Eur. Heart J. 17(3), 354
% (1996).
%
% Brennan et al., "Do existing measures of Poincare plot geometry reflect
% nonlinear features of heart rate variability?", IEEE T. Bio.-Med. Eng.
% 48(11), 1342 (2001).
%
% ---NOTES:
% Code is heavily derived from that provided by Max A. Little:
% http://www.maxlittle.net/

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

% Standard defaults
diffy = diff(y);
N = length(y); % time-series length

% ------------------------------------------------------------------------------
% Calculate pNNx percentage
% ------------------------------------------------------------------------------
% pNNx: recommendation as per Mietus et. al. 2002, "The pNNx files: ...", Heart
% strange to do this for a z-scored time series...

Dy = abs(diffy);

% Anonymous function to do the PNNx calcualtion:
% proportion of difference magnitudes greater than X*sigma
PNNxfn = @(x) mean(Dy > x / 1000);

out.pnn5  = PNNxfn(5); % 0.005*sigma
out.pnn10 = PNNxfn(10); % 0.01*sigma
out.pnn20 = PNNxfn(20); % 0.02*sigma
out.pnn30 = PNNxfn(30); % 0.03*sigma
out.pnn40 = PNNxfn(40); % 0.04*sigma

% ------------------------------------------------------------------------------
% Calculate PSD
% ------------------------------------------------------------------------------
% [Pxx, F] = psd(series,1024,1,hanning(1024),512);
[Pxx, F] = periodogram(y, hann(N)); % periodogram with hanning window

% ------------------------------------------------------------------------------
% Calculate spectral measures such as subband spectral power percentage, LF/HF ratio etc.
% ------------------------------------------------------------------------------
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
total     = fbinsize * sum(Pxx);
out.vlf   = vlfp / total * 100;
out.lf    = lfp / total * 100;
out.hf    = hfp / total * 100;

% ------------------------------------------------------------------------------
% Triangular histogram index
% ------------------------------------------------------------------------------
numBins = 10;
out.tri = length(y) / max(histcounts(y, numBins));

% ------------------------------------------------------------------------------
% Poincare plot measures:
% ------------------------------------------------------------------------------
% cf. "Do Existing Measures ... ", Brennan et. al. (2001), IEEE Trans Biomed Eng 48(11)
rmssd = std(diffy); % std of differenced series
sigma = std(y); % should be 1 for zscored time series
out.SD1 = 1 / sqrt(2) * rmssd * 1000;
out.SD2 = sqrt(2 * sigma^2 - (1 / 2) * rmssd^2) * 1000;

end
