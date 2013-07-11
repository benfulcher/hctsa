function out = WL_coeffs(y, wname, level)
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

if wmaxlev(N,wname) < level
    disp('Chosen level is too large for this wavelet on this signal');
    out = NaN; return
end

% 1. Recover a noisy signal by suppressing an
% approximation.

%% Perform a single-level wavelet decomposition 
[c,l] = wavedec(y,level,wname);

% Reconstruct detail
det = wrcoef('d',c,l,wname,level); % detail this level

det_s = sort(abs(det),'descend'); % sorted detail coefficient magnitudes

% plot(det_s);

%% Return statistics
N = length(y);
out.mean_coeff = mean(det_s);
out.max_coeff = max(det_s);
out.med_coeff = median(det_s);

% decay rate stats ('where below _ maximum' = 'wb_m')
out.wb99m = find(det_s<0.99*max(det_s),1,'first')/N;
if isempty(out.wb99m), out.wb99m = NaN; end
out.wb90m = find(det_s<0.90*max(det_s),1,'first')/N;
if isempty(out.wb90m), out.wb90m = NaN; end
out.wb75m = find(det_s<0.75*max(det_s),1,'first')/N;
if isempty(out.wb75m), out.wb75m = NaN; end
out.wb50m = find(det_s<0.50*max(det_s),1,'first')/N;
if isempty(out.wb50m), out.wb50m = NaN; end
out.wb25m = find(det_s<0.25*max(det_s),1,'first')/N;
if isempty(out.wb25m), out.wb25m = NaN; end
out.wb10m = find(det_s<0.10*max(det_s),1,'first')/N;
if isempty(out.wb10m), out.wb10m = NaN; end
out.wb1m = find(det_s<0.01*max(det_s),1,'first')/N;
if isempty(out.wb1m), out.wb1m = NaN; end

end