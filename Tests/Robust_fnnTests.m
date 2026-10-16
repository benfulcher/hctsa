classdef Robust_fnnTests < matlab.unittest.TestCase
    % Tests for NL_FNN / TISEAN false_nearest: for a scalar series the delay
    % sets the lag between embedding coordinates.

    methods (TestClassSetup)
        function setPaths(~)
            here = fileparts(fileparts(mfilename('fullpath')));
            setenv('PATH', [fullfile(here, 'Toolboxes', 'Tisean_3.0.1', 'bin'), ':', getenv('PATH')]);
        end
    end

    methods (Test)
        function delayChangesTheFalseNeighborCurve(tc)
            % x of the Lorenz system sampled every 0.02 time units: successive
            % samples are highly redundant, so a lag-1 embedding looks nearly
            % 2-dimensional, while longer lags need more dimensions (the
            % false-neighbor fraction stays high for m = 2,...,5).
            x = lorenzX(3000);
            f1 = NL_FNN(x, 1, 5, 50, 0, [], 5);
            f10 = NL_FNN(x, 10, 5, 50, 0, [], 5);
            f30 = NL_FNN(x, 30, 5, 50, 0, [], 5);
            tc.verifyLessThan(f1.pfnn_2, 0.05);
            tc.verifyLessThan(f1.pfnn_3, 0.02);
            tc.verifyGreaterThan(f10.pfnn_2, 0.1);
            tc.verifyGreaterThan(f30.pfnn_2, 0.3);
            tc.verifyGreaterThan(f30.pfnn_4, 0.1);
            tc.verifyGreaterThan(f30.pfnn_2, f10.pfnn_2);
        end

        function delayOneIsTheConsecutiveSampleEmbedding(tc)
            % -d1 is TISEAN's default, so spelling it out changes nothing, and
            % the lag-1 curve is that of the original program (values on a Lorenz
            % x series: 0.95, 0.01, 0.002, 0, 0).
            x = lorenzX(3000);
            fA = tiseanFNN(x, 1, 5, 50);
            fB = tiseanFNN(x, [], 5, 50);
            tc.verifyEqual(fA, fB, 'AbsTol', 0);
            tc.verifyEqual(fA(:,2), [0.95; 0.01; 0.002; 0; 0], 'AbsTol', 0.03);
        end

        function laggedEmbeddingEqualsLagOneOnRepeatedSeries(tc)
            % A delay-d embedding of a series in which each value is repeated d
            % times is the lag-1 embedding of the original series (the repeated
            % points coincide, and TISEAN skips zero-distance neighbors).
            z = lorenzX(1500);
            d = 4;
            x = repelem(z(:), d);
            fd = NL_FNN(x, d, 4, 0, 0, [], 5);
            f1 = NL_FNN(z, 1, 4, 0, 0, [], 5);
            tc.verifyEqual(fd.pfnn_2, f1.pfnn_2, 'AbsTol', 0.01);
            tc.verifyEqual(fd.pfnn_3, f1.pfnn_3, 'AbsTol', 0.01);
            tc.verifyEqual(fd.pfnn_1, f1.pfnn_1, 'AbsTol', 0.01);
        end

        function dimensionNeedsEnoughPointsBehindTheFraction(tc)
            % White noise: the points that find a neighbor within TISEAN's search
            % range dwindle with the dimension, and the last fractions rest on a
            % handful of points (often exactly 0). Such a fraction cannot rule a
            % dimension in or out, so the dimension is the first one with fewer
            % than ceil(3/0.05) = 60 points (a lower bound), not a "converged" one.
            rng(1, 'twister'); x = randn(1000, 1);
            M = tiseanFNN(x, 1, 10, 1); % dimension, fraction, number of points
            tc.verifyEqual(size(M, 2), 3);
            tc.verifyTrue(all(M(:,3) >= 1) && all(M(:,3) <= numel(x)));
            firstLow = find(M(:,3) < 60, 1);
            tc.verifyNotEmpty(firstLow);
            tc.verifyTrue(all(M(1:firstLow-1,2) >= 0.05));
            m = NL_FNN(x, 1, 10, 1, 1, 0.05, 5);
            tc.verifyEqual(m, M(firstLow,1));
            % 0.05 is the default threshold
            tc.verifyEqual(NL_FNN(x, 1, 10, 1, 1, [], 5), m);
            % firstunder* use the same points and never decrease as the threshold falls
            f = NL_FNN(x, 1, 10, 1, 0, [], 5);
            fu = [f.firstunder09, f.firstunder08, f.firstunder07, f.firstunder06, f.firstunder05, ...
                f.firstunder04, f.firstunder03, f.firstunder02, f.firstunder01, f.firstunder005];
            tc.verifyTrue(all(diff(fu) >= 0));
            tc.verifyEqual(f.firstunder005, m);
        end

        function lorenzDimensionUnaffectedByTheCountRule(tc)
            % A low-dimensional flow keeps many points with neighbors at every
            % dimension, so the count rule changes nothing: 3 at 0.05.
            x = lorenzX(3000);
            tc.verifyEqual(NL_FNN(x(1:5:end), 1, 10, {'ac', 1}, 1, 0.05, 5), 3);
        end
    end
end

function x = lorenzX(n)
    % Lorenz system (sigma = 10, rho = 28, beta = 8/3), fourth-order Runge-Kutta,
    % step 0.01, sampled every 2 steps; deterministic.
    s = [1; 1; 1]; h = 0.01;
    f = @(s) [10*(s(2)-s(1)); s(1)*(28-s(3))-s(2); s(1)*s(2)-8/3*s(3)];
    for i = 1:5000, s = rk4(f, s, h); end % transient
    x = zeros(n, 1);
    for i = 1:n
        s = rk4(f, s, h); s = rk4(f, s, h);
        x(i) = s(1);
    end
end

function s = rk4(f, s, h)
    k1 = f(s); k2 = f(s + h/2*k1); k3 = f(s + h/2*k2); k4 = f(s + h*k3);
    s = s + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

function M = tiseanFNN(y, d, maxm, theiler)
    % run false_nearest directly; columns: dimension, fraction false, and the
    % number of points behind the fraction
    f = BF_WriteTempFile(y);
    if isempty(d), dopt = ''; else, dopt = sprintf('-d%u ', d); end
    [~, res] = BF_TiseanSystem(sprintf('false_nearest %s-m1 -M1,%u -t%u -f5 -V0 %s', dopt, maxm, theiler, f));
    c = textscan(res, '%f%f%f%f%f');
    M = [c{1}, c{2}, c{5}];
end
