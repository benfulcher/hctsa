function out = WL_cwt(y, wname, maxScale)
% WL_cwt    Continuous wavelet transform of a time series
%
% Uses the function cwt from Matlab's Wavelet Toolbox.
%
%---INPUTS:
% y, the input time series
%
% wname, the wavelet name, e.g., 'db3' (Daubechies wavelet), 'sym2' (Symlet),
%                           etc. (see Wavelet Toolbox Documentation for all
%                           options)
%
% maxScale, the maximum scale of wavelet analysis.
%
%---OUTPUTS: statistics on the coefficients, entropy, and results of
% coefficients summed across scales.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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

% ------------------------------------------------------------------------------
% Check that a Wavelet Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('wavelet_toolbox')

% ------------------------------------------------------------------------------
%% Check inputs:
% ------------------------------------------------------------------------------
doplot = 0; % plot outputs to figures
N = length(y); % length of the time series

if nargin < 2 || isempty(wname)
    wname = 'db3';
    fprintf(1,'Using default wavelet ''%s''\n',wname);
end
if nargin < 3 || isempty(maxScale)
    maxScale = 32;
    fprintf(1,'Using default maxScale of %u\n',maxScale);
end

scales = (1:maxScale);
coeffs = cwt(y, scales, wname);

S = abs(coeffs.*coeffs); % power
SC = 100*S./sum(S(:)); % scaled power is length-dependent

% These SC values (percentage of energy in each coefficient) are what are
% displayed in a scalogram (c.f., wscalogram function)

if doplot
    figure('color','w'); box('on');
    subplot(3,1,1)
    plot(y);
    subplot(3,1,2:3);
    pcolor(SC); shading interp;
end

% ------------------------------------------------------------------------------
%% Get statistics from CWT
% ------------------------------------------------------------------------------
numEntries = size(coeffs,1)*size(coeffs,2); % number of entries in coeffs matrix

% 1) Coefficients, coeffs
allCoeffs = coeffs(:);
out.meanC = mean(allCoeffs);
out.meanabsC = mean(abs(allCoeffs));
out.medianabsC = median(abs(allCoeffs));
out.maxabsC = max(abs(allCoeffs));
out.maxonmeanC = out.maxabsC/out.meanabsC;

% 2) Power, SC -- it's highly length-dependent
% out.meanSC = mean(SC(:)); % (reproduces the mean power of power spectrum)
% out.medianSC = median(SC(:));
% out.maxSC = max(SC(:));
out.maxonmeanSC = max(SC(:))/mean(SC(:));

% Proportion of coeffs matrix over ___ maximum (thresholded)
poverfn = @(x) sum(SC(SC > x*max(SC(:))))/numEntries;
out.pover99 = poverfn(0.99);
out.pover98 = poverfn(0.88);
out.pover95 = poverfn(0.95);
out.pover90 = poverfn(0.90);
out.pover80 = poverfn(0.80);

% Distribution of scaled power
% Fit using Statistics Toolbox

if doplot
    figure('color','w');
    ksdensity(SC(:));
end

gamma_phat = gamfit(SC(:));
out.gam1 = gamma_phat(1);
out.gam2 = gamma_phat(2);

% ------------------------------------------------------------------------------
%% 2D entropy
% ------------------------------------------------------------------------------
% turn into probabilities
SC_a = SC./sum(SC(:));
% compute entropy
SC_a = SC_a(:);
out.SC_h = -sum(SC_a.*log(SC_a));

% ------------------------------------------------------------------------------
%% Weird 2-D entropy idea -- first discretize
% ------------------------------------------------------------------------------
% (i) Discretize the space into numBoxes boxes along the time axis
% Many choices, let's discretize into maximum energy
% (could also do average, or proportion inside box with more energy than
% average, ...)
numBoxes = 10;
dd_SC = zeros(maxScale, numBoxes);
cutoffs = round(linspace(0, N, numBoxes+1));
for i = 1:maxScale
   for j = 1:numBoxes
       dd_SC(i,j) = max(SC(i,cutoffs(j)+1:cutoffs(j+1)));
   end
end

% Turn into probabilities
dd_SC = dd_SC./sum(dd_SC(:));

% Compute entropy
dd_SCO = dd_SC(:);
out.dd_SC_h = -sum(dd_SCO.*log(dd_SCO));


% ------------------------------------------------------------------------------
%% Sum across scales
% ------------------------------------------------------------------------------
SSC = sum(SC);
out.max_ssc = max(SSC);
out.min_ssc = min(SSC);
out.maxonmed_ssc = max(SSC)/median(SSC);
out.pcross_maxssc50 = sum(BF_sgnchange(SSC-0.5*max(SSC))) / (N-1);
out.std_ssc = std(SSC);

% ------------------------------------------------------------------------------
%% Stationarity
% ------------------------------------------------------------------------------
% 2-way
SC_1 = SC(:,1:floor(N/2)); % collapse across scales, first half
SC_2 = SC(:,floor(N/2)+1:end); % collapse across scales, second half

mean2_1 = mean(SC_1(:));
mean2_2 = mean(SC_2(:));

std2_1 = std(SC_1(:));
std2_2 = std(SC_2(:));

% out.stat_2_m_m = mean([mean2_1 mean2_2])/mean(SC(:));
out.stat_2_m_s = mean([std2_1, std2_2])/mean(SC(:));
out.stat_2_s_m = std([mean2_1, mean2_2])/std(SC(:));
out.stat_2_s_s = std([std2_1, std2_2])/std(SC(:));

% 5-way
% I know this is terribly inefficient compared to using matrix reshape
SCs = cell(5,1);
SCs{1} = SC(:,1:floor(N/5));
SCs{2} = SC(:,floor(N/5)+1:floor(2*N/5));
SCs{3} = SC(:,floor(2*N/5)+1:floor(3*N/5));
SCs{4} = SC(:,floor(3*N/5)+1:floor(4*N/5));
SCs{5} = SC(:,floor(4*N/5)+1:end);

for i = 1:5
    out.(sprintf('mean5_%u',i)) = mean(SCs{i}(:));
    out.(sprintf('std5_%u',i)) = std(SCs{i}(:));
end

out.stat_5_m_s = mean([out.std5_1, out.std5_2, out.std5_3, out.std5_4, out.std5_5])/mean(SC(:));
out.stat_5_s_m = std([out.mean5_1, out.mean5_2, out.mean5_3, out.mean5_4, out.mean5_5])/std(SC(:));
out.stat_5_s_s = std([out.std5_1, out.std5_2, out.std5_3, out.std5_4, out.std5_5])/std(SC(:));

end
