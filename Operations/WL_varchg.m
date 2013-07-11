function out = WL_varchg(y, wname, level, maxnchpts, mindelay)
% Finds variance change points
% e.g., WL_varchg(y, 'db3', 3, 12, 0.01) uses the time series y, analyzes
% it with wavelet db3 at level 3, and finds up to a maximum of 12 change
% points with at least 1% the length of the time series separating each of
% them.
% Ben Fulcher 23/1/2010

%% Check Inputs
N = length(y);

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
if mindelay > 0 && mindelay < 1
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
    [pts_Opt, kopt, t_est] = wvarchg(det, maxnchpts, mindelay);
catch emsg
    if strcmp(emsg.identifier,'MATLAB:nomem')
       error('Not enough memory.');
    end
end

% return the number of change points
out = kopt;

% Estimated change points are close to the true change 
% points [200,600].

end