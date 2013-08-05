% SY_SlidingWindow
% 
% This function is based on sliding a window along the time series, measuring
% some quantity in each window, and outputting some summary of this set of local
% estimates of that quantity.
% 
% Another way of saying it: calculate 'windowstat' in each window, and computes
% 'acrosswindowstat' for the set of statistics calculated in each window.
% 
% INPUTS:
% 
% y, the input time series
% 
% windowstat, the measure to calculate in each window:
%               (i) 'mean', mean
%               (ii) 'std', standard deviation
%               (iii) 'ent', distribution entropy
%               (iv) 'mom3', skewness
%               (v) 'mom4', kurtosis
%               (vi) 'mom5', the fifth moment of the distribution
%               (vii) 'lillie', the p-value for a Lilliefors Gaussianity test
%               (viii) 'AC1', the lag-1 autocorrelation
%               (ix) 'apen', Approximate Entropy
% 
% acrosswindowstat, controls how the obtained sequence of local estimates is
%                   compared (as a ratio to the full time series):
%                       (i) 'std': standard deviation
%                       (ii) 'ent' histogram entropy
%                       (iii) 'apen': Approximate Entropy, ApEn(1,0.2)
%                               cf. "Approximate entropy as a measure of system
%                               complexity", S. M. Pincus, P. Natl. Acad. Sci.
%                               USA 88(6) 2297 (1991)
% 
% nseg, the number of segments to divide the time series up into, thus
%       controlling the window length
% 
% nmov, the increment to move the window at each iteration, as 1/fraction of the
%       window length (e.g., nmov = 2, means the window moves half the length of the
%       window at each increment)
% 

function out = SY_SlidingWindow(y,windowstat,acrosswindowstat,nseg,nmov)
% Ben Fulcher, 2009

doplot = 0; % plot outputs

if nargin < 2 || isempty(windowstat)
    windowstat = 'mean'; % do a sliding window mean
end
if nargin < 3 || isempty(acrosswindowstat)
    acrosswindowstat = 'std';
end
if nargin < 4 || isempty(nseg)
    nseg = 5;
end
if nargin < 5 || isempty(nmov)
    nmov = 2;
end

wlen = floor(length(y)/nseg); % size of window
inc = floor(wlen/nmov); % increment to move at each step
if inc == 0; inc = 1; end % increment rounded down to zero, prop it up

nsteps = (floor((length(y)-wlen)/inc)+1);
qs = zeros(nsteps,1);

switch windowstat
    case 'mean' % Sliding window mean
        for i = 1:nsteps
            qs(i) = mean(y((i-1)*inc + 1:(i-1)*inc + wlen));
        end
    case 'std' % Sliding window std
        for i = 1:nsteps
            qs(i) = std(y((i-1)*inc + 1:(i-1)*inc + wlen));
        end
    case 'ent' % Sliding window distributional entropy
        for i = 1:nsteps
            ksstats = DN_FitKernelSmooth(y((i-1)*inc + 1:(i-1)*inc + wlen),'entropy');
            qs(i) = ksstats.entropy;
        end
    case 'apen' % Sliding window ApEn
        for i = 1:nsteps
            qs(i) = EN_ApEn(y((i-1)*inc + 1:(i-1)*inc + wlen),1,0.2);
        end
    case 'mom3' % Third moment
        for i = 1:nsteps
            qs(i) = DN_Moments(y((i-1)*inc + 1:(i-1)*inc + wlen),3);
        end
    case 'mom4' % Fourth moment
        for i = 1:nsteps
            qs(i) = DN_Moments(y((i-1)*inc + 1:(i-1)*inc + wlen),4);
        end
    case 'mom5' % Fifth moment
        for i = 1:nsteps
            qs(i) = DN_Moments(y((i-1)*inc + 1:(i-1)*inc + wlen),5);
        end
    case 'lillie' % Lilliefors test
        for i = 1:nsteps
            qs(i) = HT_DistributionTest(y((i-1)*inc + 1:(i-1)*inc + wlen),'lillie','norm');
        end
    case 'AC1' % Lag-1 autocorrelation
        for i = 1:nsteps
            qs(i) = CO_AutoCorr(y((i-1)*inc + 1:(i-1)*inc + wlen),1);
        end
    otherwise
        error('Unknown statistic ''%s''',windowstat)
end

if doplot
    figure('color','w'); box('on');
    plot(round(wlen/2):inc:(nsteps-1)*inc+round(wlen/2),qs,'r');
end

switch acrosswindowstat
    case 'std'
        out = std(qs)/std(y);
    case 'apen'
        out = EN_ApEn(qs,1,0.2); % ApEn of the sliding window measures
    case 'ent'
        kssimpouts = DN_FitKernelSmooth(qs); % get a load of statistics from kernel-smoothed distribution
        out = kssimpouts.entropy; % distributional entropy
    otherwise
        error('Unknown statistic: ''%s''',acrosswindowstat)
end

end