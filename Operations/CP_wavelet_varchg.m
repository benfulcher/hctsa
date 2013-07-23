% CP_wavelet_varchg
% 
% Finds variance change points using functions from Matlab's Wavelet Toolbox,
% including the primary function wvarchg, which estimates the change points in
% the time series.
% 
% INPUTS:
% 
% y, the input time series
% 
% wname, the name of the mother wavelet to analyze the data with: e.g., 'db3',
%           'sym2', cf. Wavelet Toolbox Documentation for details
% 
% level, the level of wavelet decomposition
% 
% maxnchpts, the maximum number of change points
% 
% mindelay, the minimum delay between consecutive change points (can be
%           specified as a proportion of the time-series length, e.g., 0.02
%           ensures that change points are separated by at least 2% of the
%           time-series length)
% 
% The output from this function is the optimal number of change points.
% 

function out = CP_wavelet_varchg(y, wname, level, maxnchpts, mindelay)
% Ben Fulcher, 23/1/2010

%% Check Inputs
N = length(y); % time-series length

if nargin < 2 || isempty(wname)
    wname = 'db3'; % default wavelet
end

if nargin < 3 || isempty(level)
   level = 3; % level of wavelet decomposition
end
if strcmp(level,'max')
    level = wmaxlev(N,wname);
end

if nargin < 4 || isempty(maxnchpts)
   maxnchpts = 5; % maximum number of change points
end

if nargin < 5 || isempty(mindelay)
     mindelay = 0.01; % 1% of the time series length
end
if (mindelay > 0) && (mindelay < 1)
   mindelay = ceil(mindelay*N);
end

if wmaxlev(N, wname) < level
    error('Chosen level, %u, is too large for this wavelet on this signal. Sorry.', level);
end

% The aim of this example is to recover the
% change points in signal y.
% In addition, this example illustrates how the GUI
% tools propose change point locations for interval
% dependent de-noising thresholds.
% 1. Recover a noisy signal by suppressing an
% approximation.

%% Perform a single-level wavelet decomposition 
[c, l] = wavedec(y,level,wname);

% Reconstruct detail at the same level.
det = wrcoef('d',c,l,wname,level);

% % 2. Replace 2% of the greatest (absolute) values by the mean
% % in order to remove almost all the signal.
x = sort(abs(det));
v2p100 = x(fix(length(x)*0.98));
det(abs(det)>v2p100) = mean(det);

% keyboard

% 3. Use wvarchg for estimating the change points
try
    [~, kopt, ~] = wvarchg(det, maxnchpts, mindelay);
catch emsg
    if strcmp(emsg.identifier,'MATLAB:nomem')
       error('Not enough memory.');
    end
end

% return the number of change points
out = kopt;

end