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
[c, l] = wavedec(y,level,wname);

% Reconstruct detail
det = wrcoef('d',c,l,wname,level); % detail this level

det_s = sort(abs(det),'descend'); % sorted detail coefficient magnitudes

% plot(det_s);

%% Return statistics
out.mean_coeff = mean(det_s);
out.max_coeff = max(det_s);
out.med_coeff = median(det_s);

% decay rate stats ('where below _ maximum' = 'wb_m')
out.wb99m = findmythreshold(0.99);
out.wb90m = findmythreshold(0.90);
out.wb75m = findmythreshold(0.75);
out.wb50m = findmythreshold(0.50);
out.wb25m = findmythreshold(0.25);
out.wb10m = findmythreshold(0.10);
out.wb1m = findmythreshold(0.01);

    function i = findmythreshold(x)
        % where drops below proportion x of maximum
        i = find(det_s < x*max(det_s),1,'first')/N;
        if isempty(i), i = NaN; end
    end

end