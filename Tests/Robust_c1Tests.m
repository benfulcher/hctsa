classdef Robust_c1Tests < matlab.unittest.TestCase
    % Tests of NL_c1 (TISEAN c1/c2d): the patched c1 must terminate where stock TISEAN
    % searches for neighbors for ever. These need the compiled TISEAN binaries
    % (Toolboxes/compile_tisean.m); with an unpatched c1 they would block until the
    % 600 s timeout of BF_TiseanSystem.

    methods (Test)
        function powerOfTwoNumberOfEmbeddedPointsTerminates(tc)
            % 520 samples with delay 2 give 512 embedded points at m = 5: the fixed mass reaches 1
            rng(1); y = zscore(cumsum(randn(520, 1)));
            t = tic;
            o = NL_c1(y, 2, [1, 7], 26, 156);
            tc.verifyLessThan(toc(t), 120);
            tc.verifyTrue(isstruct(o));
            tc.verifyTrue(isfinite(o.bestestd));
        end

        function lengthsNearMultipleOf128AreNotTrimmed(tc)
            % the former workaround dropped the last mod(N,128)+1 samples; the result must now
            % use all of them, i.e. differ from the result for the trimmed series
            rng(2); y = zscore(randn(640, 1));
            a = NL_c1(y, 1, [2, 5], 0.02, 0.5);
            b = NL_c1(y(1:end-1), 1, [2, 5], 0.02, 0.5);
            tc.verifyTrue(isstruct(a) && isstruct(b));
            tc.verifyNotEqual(a.bestestd, b.bestestd);
        end

        function allPointsAsCentersTerminates(tc)
            % more centers requested (Nref = N) than there are embedded points
            rng(3); y = zscore(cumsum(randn(300, 1)));
            o = NL_c1(y, 2, [2, 6], 0.02, 1);
            tc.verifyTrue(isstruct(o) || isnan(o));
        end

        function timeSeparationWithoutNeighborsIsNaN(tc)
            rng(4); y = randn(300, 1);
            tc.verifyTrue(isnan(NL_c1(y, 1, [2, 4], 200, 0.5)));
        end

        function dimensionIndependentOfNumberOfReferencePoints(tc)
            % The mean log radius is averaged over the reference points actually used. It used to be
            % divided by (Nref - (m-1)*delay), which with few reference points and a long delay
            % inflated it and pulled the dimension estimate down (about 30% at m = 4 for 100
            % reference points of iid noise, delay 10). Estimates must agree across Nref.
            rng(7); y = randn(2000, 1);
            Nrefs = [100, 300, 2000];
            maxmd = zeros(size(Nrefs));
            for i = 1:numel(Nrefs)
                o = NL_c1(y, 10, [2, 4], 0.02, Nrefs(i));
                maxmd(i) = o.maxmd; % the m = 4 estimate (true value 4 for iid noise)
            end
            tc.verifyLessThan(max(maxmd)/min(maxmd), 1.07);
            tc.verifyGreaterThan(min(maxmd), 3.7);
        end

        function deterministic(tc)
            rng(5); y = zscore(cumsum(randn(1000, 1)));
            a = NL_c1(y, 1, [1, 7], 0.02, 0.5);
            b = NL_c1(y, 1, [1, 7], 0.02, 0.5);
            tc.verifyEqual(a, b);
        end
    end
end
