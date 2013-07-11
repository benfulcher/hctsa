function out = WL_scal2frq(y, wname, amax, delta)
% Ben Fulcher 26/1/2010. Adapted from example in MATLAB wavelet toolbox
% It's a bit stupid really.

%% Check Inputs
N = length(y); % length of time series

if nargin < 2 || isempty(wname)
    disp('Wavelet not specified -- using the default db3 wavelet')
    wname = 'db3';
end

if nargin < 3 || isempty(amax)
    amax = 5; % maximum 'scale'
end
maxlevel = wmaxlev(N,wname); % maximum level for this time-series length
if strcmp(amax,'max') % set to maximum for this wavelet
    amax = wmaxlev(N,wname);
end

if nargin < 4 || isempty(delta)
    delta = 1; % the sampling period
end

if maxlevel < amax
    fprintf(1,'Chosen level (%u) is too large for this wavelet on this signal...',amax);
    amax = maxlevel;
    fprintf(1,' changed to maximum level computed with wmaxlev: %u\n',amax);
end

%% Do your thing.
% This example demonstrates that, starting from the periodic function
% x(t) = 5*sin(5t) + 3*sin(2t) + 2*sin(t), the scal2frq function translates
% the scales corresponding to the maximum values of the CWT coefficients
% to pseudo-frequencies ([0.796 0.318 0.159]), which are near to the true
% frequencies ([5 2 1] / (2*pi) =~ [0.796 0.318 0.159]).

% delta = 0.1;
% wname = 'coif3';

% Define scales.
scales = 1:amax;
a = 2.^scales;

% Compute associated pseudo-frequencies.
f = scal2frq(a, wname, delta);

% Compute associated pseudo-periods.
per = 1./f; 

% Decompose the time series at level specified as maximum
[c, l] = wavedec(y, amax, wname);

% Estimate standard deviation of detail coefficients.
stdc = wnoisest(c,l,scales);

% plot(stdc) % plot them

% Compute identified period.
[~, jmax] = max(stdc); % level with highest energy coefficients

out.lmax = jmax; % level with highest energy coefficients
out.period = per(jmax); % output dominant period
out.pf = f(jmax); % output dominant pseudo-frequency


      
end