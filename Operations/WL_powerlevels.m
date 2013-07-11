function out = WL_powerlevels(y, wname, maxlevel)
% Compares the detail coefficients at each level
% Ben Fulcher 23/1/2010

doplot = 0; % can plot outputs

%% Check Inputs
N = length(y);

if nargin < 2 || isempty(wname)
    wname = 'db3'; % default wavelet
end
if nargin < 3 || isempty(maxlevel)
   maxlevel = 20; % maximum wavelet decomposition level
end
if strcmp(maxlevel,'max')
    maxlevel = wmaxlev(N,wname);
end

if wmaxlev(N, wname) < maxlevel
    fprintf(1,'Chosen wavelet level is too large for the %s wavelet for this signal of length N = %u\n',wname,N);
    maxlevel = wmaxlev(N,wname);
    fprintf(1,'Using a wavelet level of %u instead.\n',maxlevel)
end


%% Perform a single-level wavelet decomposition
means = zeros(maxlevel,1); % mean detail coefficient magnitude at each level
medians = zeros(maxlevel,1); % median detail coefficient magnitude at each level
maxs = zeros(maxlevel,1); % max detail coefficient magnitude at each level

for k = 1:maxlevel
    level = k;
    
    [c, l] = wavedec(y,level,wname);
    % Reconstruct detail at this level
    det = wrcoef('d',c,l,wname,level);
    
    means(k) = mean(abs(det));
    medians(k) = median(abs(det));
    maxs(k) = max(abs(det));
end

%% Plot the bad boy

if doplot
    subplot(5,1,1:2); title('signal')
    plot(y);
    subplot(5,1,3); title('means');
    plot(means)
    subplot(5,1,4); title('medians');
    plot(medians)
    subplot(5,1,5); title('maxs');
    plot(maxs);
end

%% Return statistics on detail coefficients
% Sort
means_s = sort(means,'descend');
medians_s = sort(medians,'descend');
maxs_s = sort(maxs,'descend');

% What is the maximum across these levels
out.max_mean = means_s(1);
out.max_median = medians_s(1);
out.max_max = maxs_s(1);

% stds
out.std_mean = std(means);
out.std_median = std(medians);
out.std_max = std(maxs);

% At what level is the maximum
out.wheremax_mean = find(means == means_s(1));
out.wheremax_median = find(medians == medians_s(1));
out.wheremax_max = find(maxs == maxs_s(1));

% Size of maximum (relative to next maximum)
out.max1on2_mean = means_s(1)/means_s(2);
out.max1on2_median = medians_s(1)/medians_s(2);
out.max1on2_max = maxs_s(1)/maxs_s(2);

% Where sum of values to left equals sum of values to right
% Measure of centrality
out.wslesr_mean = SUB_slosr(means);
out.wslesr_median = SUB_slosr(medians);
out.wslesr_max = SUB_slosr(maxs);

% What's the correlation between maximum and median
r = corrcoef(maxs,medians);
out.corrcoef_max_medians = r(1,2);


function meisgorilla = SUB_slosr(xx)
    maxlevel = length(xx);
    slosr = zeros(maxlevel-2,1);
    for i = 2:maxlevel-1
        slosr(i-1) = sum(xx(1:i-1))/sum(xx(i+1:end));
    end
    absm1 = abs(slosr-1); % how close to 1 (the same sum on either side) each is
    meisgorilla = find(absm1 == min(absm1),1,'first') + 1;
end
    
end