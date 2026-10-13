classdef Robust_sysidTests < matlab.unittest.TestCase
    % Tests of hctsa's own linear system-identification helpers (BF_SubspaceID, BF_SSPredict,
    % BF_ARMAfit, BF_ARfb, BF_ARorderAIC, BF_Pinv) and of the MF_ operations built on them.

    methods (Test)
        function diffusePredictionsOfAnARModel(tc)
            % AR(1): the first sample is not predicted; afterwards yhat(t) = 0.6 y(t-1)
            rng(1); y = filter(1, [1, -0.6], randn(300, 1));
            [A, K, C] = BF_ARMAtoSS([1, -0.6], 1);
            yp = BF_SSPredict(A, K, C, y, 1);
            tc.verifyTrue(isnan(yp(1)));
            tc.verifyEqual(yp(2:end), 0.6 * y(1:end-1), 'AbsTol', 1e-12);
            yp2 = BF_SSPredict(A, K, C, y, 2); % two steps ahead: 0.36 y(t-2)
            tc.verifyTrue(all(isnan(yp2(1:2))));
            tc.verifyEqual(yp2(3:end), 0.36 * y(1:end-2), 'AbsTol', 1e-12);
        end

        function diffusePredictionsIgnoreTheLevel(tc)
            % a random walk with a large offset is predicted by its previous value from sample 2
            rng(2); y = 1e6 + cumsum(randn(200, 1));
            yp = BF_SSPredict(1, 1, 1, y, 1);
            tc.verifyEqual(yp(2:end), y(1:end-1), 'RelTol', 1e-12);
        end

        function subspaceRecoversAnAR2(tc)
            rng(3);
            r = 0.9; w = 0.5;
            y = filter(1, [1, -2*r*cos(w), r^2], randn(8000, 1));
            m = BF_SubspaceID(zscore(y), 2);
            lam = sort(abs(angle(eig(m.A))));
            tc.verifyEqual(abs(eig(m.A)), [r; r], 'AbsTol', 0.03);
            tc.verifyEqual(lam, [w; w], 'AbsTol', 0.03);
            tc.verifyTrue(all(m.C >= 0)); % the sign convention
            % the loss is the innovation variance (of the z-scored series)
            tc.verifyEqual(m.loss, 1 / var(y), 'RelTol', 0.05);
        end

        function subspaceBestOrderOfAnExactSinusoid(tc)
            % the AIC (with its rounding-level floor) picks the lowest order that is exact: a
            % (zero-mean) sinusoid is an exact AR(2) process / an order-2 state-space model
            t = (1:600)';
            y = sin(0.2 * t);
            tc.verifyEqual(BF_ARorderAIC(y), 2);
            m = BF_SubspaceID(y, 'best');
            tc.verifyEqual(m.order, 2);
        end

        function armaFitMinimizesTheConditionalSumOfSquares(tc)
            rng(4);
            y = filter([1, 0.5], [1, -0.6], randn(20000, 1));
            [a, c, sd, loss] = BF_ARMAfit(y, 1, 1);
            tc.verifyEqual(a(2), -0.6, 'AbsTol', 0.03);
            tc.verifyEqual(c(2), 0.5, 'AbsTol', 0.03);
            tc.verifyLessThan(max(sd), 0.03);
            tc.verifyEqual(loss, 1, 'AbsTol', 0.03);
            % a stationary point of the conditional sum of squares (central differences)
            y = y(1:600);
            [a, c] = BF_ARMAfit(y, 1, 1);
            css = @(th) sum(SUB_condResid(y, th(1), th(2)).^2);
            th = [a(2), c(2)]; h = 1e-6;
            g = [css(th + [h, 0]) - css(th - [h, 0]), css(th + [0, h]) - css(th - [0, h])] / (2*h);
            tc.verifyLessThan(max(abs(g)), 1e-4 * css(th));
        end

        function armaFitIsInvertible(tc)
            rng(5);
            y = filter([1, -0.9], 1, randn(4000, 1)); % zero at 0.9
            [~, c] = BF_ARMAfit(y, 1, 1);
            tc.verifyTrue(all(abs(roots(c)) < 1));
        end

        function forwardBackwardARIsTimeSymmetric(tc)
            rng(6); y = filter(1, [1, -0.6, 0.2], randn(200, 1));
            tc.verifyEqual(BF_ARfb(y, 3), BF_ARfb(flipud(y), 3), 'AbsTol', 1e-10);
        end

        function pinvCutoff(tc)
            M = [1, 0; 0, 1e-12; 0, 0];
            tc.verifyEqual(BF_Pinv(M), [1, 0, 0; 0, 0, 0]); % the 1e-12 direction is below the cutoff
            rng(7); X = randn(30, 4);
            tc.verifyEqual(BF_Pinv(X), pinv(X), 'AbsTol', 1e-12);
        end

        function exactlyPredictableResidualsGiveNaNAutocorrelation(tc)
            t = (1:500)';
            o = MF_steps_ahead(sin(0.2 * t), 'ar', 2, 3); % an exact AR(2) process
            tc.verifyTrue(isnan(o.ac1_h1));
            tc.verifyLessThan(o.stde_h1, 1e-9);
        end

        function residualSummaryOfRoundingLevelResiduals(tc)
            % residuals that do not vary beyond rounding: only meane, meanabs, stde are defined
            rng(9); y = randn(300, 1);
            for e = {1e-13 * randn(300, 1), 0.3 * ones(300, 1), zeros(300, 1)}
                o = MF_ResidualAnalysis(e{1}, y, 'full');
                tc.verifyTrue(all(isfinite([o.meane, o.meanabs, o.stde])));
                tc.verifyTrue(all(isnan([o.maxonstd, o.ac1, o.ac2, o.ac3, o.propbth, o.taurat, ...
                                o.sws, o.swm, o.ftbth, o.normksstat, o.popt, o.minsbc])));
            end
            o = MF_ResidualAnalysis(randn(300, 1), y, 'full'); % ordinary residuals: all defined
            tc.verifyTrue(all(isfinite(cell2mat(struct2cell(o)))));
        end

        function tooShortSeriesGiveNaN(tc)
            rng(8); y = randn(12, 1);
            tc.verifyTrue(isequaln(MF_StateSpace_n4sid(y, 3, 0.5, 1), NaN));
            tc.verifyTrue(isequaln(MF_armax(y, [3, 1], 0.5, 1), NaN));
        end
    end
end

function e = SUB_condResid(y, a1, c1)
    % Conditional ARMA(1,1) residuals e(t), t >= 2, with e(1) = 0
    e = zeros(size(y));
    for t = 2:length(y)
        e(t) = y(t) + a1 * y(t-1) - c1 * e(t-1);
    end
    e = e(2:end);
end
