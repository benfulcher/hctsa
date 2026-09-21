% hctsa_smoke_report
% Diagnostic (not a test): evaluates the full hctsa feature library on a
% small panel of diverse synthetic time series and reports
%   - operations that error on any series,
%   - features that are NaN/special on every series,
%   - features that take a single value across all series (candidates for
%     deregistration, to be confirmed on a real dataset such as
%     Empirical1000 before acting),
%   - the slowest master operations.
% Run after a batch of changes to Operations/ or the INP_*.txt files.
% Takes ~10 minutes (13 series x ~30-60 s each). The pass/fail robustness
% checks (no errors, determinism) live in RobustnessTests.m.

run(fullfile(fileparts(mfilename('fullpath')), '..', 'startup.m'));

% ------------------------------------------------------------------------------
% The panel: the adversarial series from RobustnessTests plus ordinary ones
% ------------------------------------------------------------------------------
names = {'quantized', 'nearConstant', 'randomWalk', 'spikes', 'short', 'ar1'};
series = cellfun(@(n) RobustnessTests.adversarialSeries(n), names, 'UniformOutput', false);
rng(20260920);
N = 1000;
series{end+1} = randn(N, 1);                                        names{end+1} = 'whiteNoise';
series{end+1} = sin(2*pi*(1:N)'/25) + 0.3*randn(N, 1);              names{end+1} = 'noisySine';
x = zeros(N, 1); x(1) = 0.3; for t = 2:N, x(t) = 4*x(t-1)*(1 - x(t-1)); end
series{end+1} = x;                                                  names{end+1} = 'logisticMap';
series{end+1} = trnd(3, N, 1);                                      names{end+1} = 'heavyTailed';
series{end+1} = (1:N)'/N*5 + randn(N, 1);                           names{end+1} = 'trend';
e = randn(N+1, 1); series{end+1} = e(2:end) - 0.9*e(1:end-1);       names{end+1} = 'ma1';
series{end+1} = mod((1:N)', 7) + 0.05*randn(N, 1);                  names{end+1} = 'sawtooth';
numSeries = numel(series);

% ------------------------------------------------------------------------------
% Compute
% ------------------------------------------------------------------------------
Operations = TS_ReadInputFile('ops', 'INP_ops_hctsa.txt', false);
MasterOperations = TS_ReadInputFile('mops', 'INP_mops_hctsa.txt', false);
[Operations, MasterOperations] = TS_LinkOperationsWithMasters(Operations, MasterOperations);
numOps = height(Operations);
D = nan(numSeries, numOps); Q = nan(numSeries, numOps); T = nan(numSeries, numOps);
for i = 1:numSeries
    fprintf(1, '[%u/%u] %s (N = %u)...', i, numSeries, names{i}, length(series{i}));
    tic
    [D(i,:), T(i,:), Q(i,:)] = TS_CalculateFeatureVector(series{i}, false, Operations, MasterOperations, true, 'fast');
    fprintf(1, ' %.1f s, %u errors.\n', toc, sum(Q(i,:) == 1));
end
codeOf = @(mid) MasterOperations.Code{MasterOperations.ID == mid};

% ------------------------------------------------------------------------------
% Report
% ------------------------------------------------------------------------------
fprintf(1, '\n=== Errors (quality = 1), by master operation ===\n');
errMasters = unique(Operations.MasterID(any(Q == 1, 1)));
for i = 1:numel(errMasters)
    onSeries = names(any(Q(:, Operations.MasterID == errMasters(i)) == 1, 2));
    fprintf(1, '  %s  <- %s\n', codeOf(errMasters(i)), strjoin(onSeries, ', '));
end
if isempty(errMasters), fprintf(1, '  (none)\n'); end

fprintf(1, '\n=== Features never returning a real value (NaN/special on every series) ===\n');
neverGood = find(all(Q > 0, 1));
[m, ~, ic] = unique(Operations.MasterID(neverGood));
for i = 1:numel(m)
    fprintf(1, '  [%u fields] %s\n', sum(ic == i), codeOf(m(i)));
end
if isempty(m), fprintf(1, '  (none)\n'); end

fprintf(1, '\n=== Features constant across all series (>= %u good values) -- deregistration candidates ===\n', numSeries - 3);
isConst = false(1, numOps);
for j = 1:numOps
    v = D(Q(:, j) == 0, j);
    isConst(j) = numel(v) >= numSeries - 3 && max(v) == min(v);
end
[m, ~, ic] = unique(Operations.MasterID(isConst));
constIdx = find(isConst);
for i = 1:numel(m)
    jj = constIdx(ic == i);
    vals = arrayfun(@(j) D(find(Q(:, j) == 0, 1), j), jj);
    fprintf(1, '  %s : %s\n', codeOf(m(i)), strjoin(compose('%s = %g', string(Operations.Name(jj)), vals'), ', '));
end
if isempty(m), fprintf(1, '  (none)\n'); end

fprintf(1, '\n=== Slowest master operations (mean seconds over the N = %u series) ===\n', N);
isN = cellfun(@length, series) == N;
mt = zeros(height(MasterOperations), 1);
for i = 1:height(MasterOperations)
    col = find(Operations.MasterID == MasterOperations.ID(i), 1);
    if ~isempty(col), mt(i) = mean(T(isN, col), 'omitnan'); end
end
[~, ord] = sort(mt, 'descend');
for i = 1:20
    fprintf(1, '  %6.2f s (%4.1f%%)  %s\n', mt(ord(i)), 100*mt(ord(i))/sum(mt), MasterOperations.Code{ord(i)});
end
fprintf(1, '  (total %.1f s per series)\n', sum(mt));
