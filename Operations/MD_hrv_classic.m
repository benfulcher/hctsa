function out = MD_hrv_classic(y)
% MD_hrv_classic   Classic heart rate variability (HRV) statistics.
%
% Typically assumes an NN/RR time series in units of seconds. Returns pNNx-style
% measures (the proportion of successive differences larger than a multiple of
% their robust standard deviation, so that they do not depend on the units or
% scale of the series), the
% proportions of power in the very-low, low and high frequency bands of a
% Hann-windowed periodogram and the ratio of low to high, the triangular
% histogram index, and Poincare plot measures. The frequency bands are applied
% to the periodogram's normalized frequency in radians per sample, not in hertz.
%
% ---INPUTS:
% y, the input time series
%
% ---OUTPUTS:
% pnnrel025, pnnrel05, pnnrel1, pnnrel2, pnnrel3: proportion of successive
%       differences of y whose magnitude exceeds 0.25, 0.5, 1, 2 and 3 times the
%       robust standard deviation of the successive differences,
%       sigD = median(|d - median(d)|)/0.6745 where d = diff(y) (0.6745 makes it
%       consistent with the standard deviation for Gaussian increments). If
%       sigD = 0 (more than half the increments equal), the mean absolute
%       deviation of the increments about their median, times sqrt(pi/2), is
%       used instead; if that is also 0 (all increments equal), NaN is returned
%       for all five.
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
%       sqrt(2*std(y)^2 - std(diff(y))^2/2), times 1000 (not registered as a
%       feature: for a z-scored series it is 1000*sqrt(1 + AC1) up to end effects, and so
%       is redundant with the lag-1 autocorrelation)
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
% The original pnn5 to pnn40 used fixed thresholds of x/1000 on the z-scored series
% (0.005 to 0.04 standard deviations), so for most series they were close to 1 and
% almost constant across series. The pnnrel* measures replace them (pnn5 -> pnnrel025,
% pnn10 -> pnnrel05, pnn20 -> pnnrel1, pnn30 -> pnnrel2, pnn40 -> pnnrel3).
%
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
% Calculate pNNx: proportion of |successive differences| exceeding c robust SDs
% ------------------------------------------------------------------------------
% pNNx: cf. Mietus et. al. 2002, "The pNNx files: ...", Heart. The fixed thresholds
% x/1000 (in the units of the series) are replaced by multiples of the robust
% standard deviation of the increments, so the measure does not depend on the units.
Dy = abs(diffy);
sigD = median(abs(diffy - median(diffy))) / 0.6745; % robust (MAD-based) SD of increments
tolD = 1e-10 * std(y); % treat a spread below this (rounding error) as zero
if sigD <= tolD
	% Over half the increments are equal (e.g., a quantized series): fall back to the
	% mean absolute deviation about the median (consistent with the SD for Gaussian
	% increments). Zero only if all the increments are equal, in which case NaN.
	sigD = mean(abs(diffy - median(diffy))) * sqrt(pi / 2);
end
if sigD <= tolD
	sigD = NaN;
end
if isnan(sigD)
	PNNxfn = @(c) NaN; % (a comparison with NaN would otherwise give 0)
else
	PNNxfn = @(c) mean(Dy > c * sigD);
end

out.pnnrel025 = PNNxfn(0.25);
out.pnnrel05 = PNNxfn(0.5);
out.pnnrel1 = PNNxfn(1);
out.pnnrel2 = PNNxfn(2);
out.pnnrel3 = PNNxfn(3);

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
out.tri = length(y) / max(histcounts(y, BF_HistEdges(y, numBins))); % equal-width bins spanning the data (BF_HistEdges)

% ------------------------------------------------------------------------------
% Poincare plot measures:
% ------------------------------------------------------------------------------
% cf. "Do Existing Measures ... ", Brennan et. al. (2001), IEEE Trans Biomed Eng 48(11)
rmssd = std(diffy); % std of differenced series
sigma = std(y); % should be 1 for zscored time series
out.SD1 = 1 / sqrt(2) * rmssd * 1000;
out.SD2 = sqrt(2 * sigma^2 - (1 / 2) * rmssd^2) * 1000;

end
