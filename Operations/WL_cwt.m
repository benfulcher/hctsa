function out = WL_cwt(y, wname, maxscale)
% Uses Wavelet, Statistics Toolboxes in MATLAB
% Ben Fulcher 26/1/2010

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

% subplot(3,1,1)
% plot(y);
% subplot(3,1,2:3);
% pcolor(SC); shading interp;
% keyboard

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

% proportion of coeffs matrix over ___ maximum (thresholded)
out.pover99 = sum(SC(SC > 0.99*max(SC(:))))/Nentries;
out.pover98 = sum(SC(SC > 0.98*max(SC(:))))/Nentries;
out.pover95 = sum(SC(SC > 0.95*max(SC(:))))/Nentries;
out.pover90 = sum(SC(SC > 0.90*max(SC(:))))/Nentries;
out.pover80 = sum(SC(SC > 0.80*max(SC(:))))/Nentries;

% Distribution of scaled power
% Fit using Statistics Toolbox

% figure('color','w');
% ksdensity(SC(:));
    
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
out.pcross_maxssc50 = length(sgnchange(SSC-0.5*max(SSC)))/N;
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