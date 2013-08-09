% WL_cwt
% 
% Applies a continuous wavelet transform to the time series using the function
% cwt from Matlab's Wavelet Toolbox.
% 
% INPUTS:
% y, the input time series
% 
% wname, the wavelet name, e.g., 'db3' (Daubechies wavelet), 'sym2' (Symlet),
%                           etc. (see Wavelet Toolbox Documentation for all
%                           options)
% 
% maxscale, the maximum scale of wavelet analysis.
% 
% 
% Outputs from this function are statistics on the coefficients, entropy, and
% results of coefficients summed across scales.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = WL_cwt(y, wname, maxscale)
% Ben Fulcher, 26/1/2010

%% Check that a Wavelet Toolbox license exists:
a = license('test','wavelet_toolbox');
if a == 0
    error('This function requires Matlab''s Wavelet Toolbox');
end
% Try to check out a license:
[lic_free,~] = license('checkout','wavelet_toolbox');
if lic_free == 0
    error('Could not obtain a license for Matlab''s Wavelet Toolbox');
end

% Check inputs:
doplot = 0; % plot outputs to figures
N = length(y); % length of the time series

if nargin < 2 || isempty(wname)
    wname = 'db3';
    fprintf(1,'Using default wavelet ''%s''\n',wname)
end
if nargin < 3 || isempty(maxscale)
    maxscale = 32;
    fprintf(1,'Using default maxscale of %u\n',maxscale)
end

scales = (1:maxscale);
coeffs = cwt(y, scales, wname);

S = abs(coeffs.*coeffs); % power
SC = 100*S./sum(S(:)); % scaled power

% These SC values (percentage of energy in each coefficient) are what are
% displayed in a scalogram (c.f., wscalogram function)

if doplot
    figure('color','w'); box('on');
    subplot(3,1,1)
    plot(y);
    subplot(3,1,2:3);
    pcolor(SC); shading interp;
end

%% Get statistics from CWT
Nentries = size(coeffs,1)*size(coeffs,2); % number of entries in coeffs matrix

% 1) Coefficients, coeffs
out.meanC = mean(coeffs(:));
out.meanabsC = mean(abs(coeffs(:)));
out.medianabsC = median(abs(coeffs(:)));
out.maxabsC = max(abs(coeffs(:)));
out.maxonmedC = max(abs(coeffs(:)))/median(abs(coeffs(:)));


% 2) Power, SC
out.meanSC = mean(SC(:));
out.medianSC = median(SC(:));
out.maxSC = max(SC(:));
out.maxonmedSC = max(SC(:))/median(SC(:));

% Proportion of coeffs matrix over ___ maximum (thresholded)
poverfn = @(x) sum(SC(SC > x*max(SC(:))))/Nentries;
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
    
out.exp_muhat = expfit(SC(:));
gamma_phat = gamfit(SC(:));
out.gam1 = gamma_phat(1);
out.gam2 = gamma_phat(2);

%% 2D entropy
% turn into probabilities
SC_a = SC./sum(SC(:));
% compute entropy
SC_a = SC_a(:);
out.SC_h = -sum(SC_a.*log(SC_a));

%% Weird 2-D entropy idea -- first discretize
% (i) Discretize the space into nboxes boxes along the time axis
% Many choices, let's discretize into maximum energy
% (could also do average, or proportion inside box with more energy than
% average, ...)
nboxes = 10;
dd_SC = zeros(maxscale, nboxes);
cutoffs = round(linspace(0, N, nboxes+1));
for i = 1:maxscale
   for j = 1:nboxes
       dd_SC(i,j) = max(SC(i,cutoffs(j)+1:cutoffs(j+1)));
   end
end

% turn into probabilities
dd_SC = dd_SC./sum(dd_SC(:));

% compute entropy
dd_SCO = dd_SC(:);
out.dd_SC_h = -sum(dd_SCO.*log(dd_SCO));


%% Sum across scales
SSC = sum(SC);
out.max_ssc = max(SSC);
out.min_ssc = min(SSC);
out.maxonmed_ssc = max(SSC)/median(SSC);
out.pcross_maxssc50 = sum(BF_sgnchange(SSC-0.5*max(SSC))) / (N-1);
out.std_ssc = std(SSC);

%% Stationarity
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
SC_1 = SC(:,1:floor(N/5));
SC_2 = SC(:,floor(N/5)+1:floor(2*N/5));
SC_3 = SC(:,floor(2*N/5)+1:floor(3*N/5));
SC_4 = SC(:,floor(3*N/5)+1:floor(4*N/5));
SC_5 = SC(:,floor(4*N/5)+1:end);

for i = 1:5
    eval(sprintf('mean5_%u = mean(SC_%u(:));',i,i))
    eval(sprintf('std5_%u = std(SC_%u(:));',i,i))
end

% mean5_1 = mean(SC_1(:));
% mean5_2 = mean(SC_2(:));
% mean5_3 = mean(SC_3(:));
% mean5_4 = mean(SC_4(:));
% mean5_5 = mean(SC_5(:));

% std5_1 = std(SC_1(:));
% std5_2 = std(SC_2(:));
% std5_3 = std(SC_3(:));
% std5_4 = std(SC_4(:));
% std5_5 = std(SC_5(:));

% out.stat_5_m_m = mean([mean5_1 mean5_2 mean5_3 mean5_4
% mean5_5])/mean(SC(:));
out.stat_5_m_s = mean([std5_1, std5_2, std5_3, std5_4, std5_5])/mean(SC(:));
out.stat_5_s_m = std([mean5_1, mean5_2, mean5_3, mean5_4, mean5_5])/std(SC(:));
out.stat_5_s_s = std([std5_1, std5_2, std5_3, std5_4, std5_5])/std(SC(:));

end