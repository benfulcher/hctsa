classdef SamplingRateGuardsTests < matlab.unittest.TestCase
    % Guards on the sampling-rate normalization: the length floor of the
    % 'decimate_ac1e' preprocessing (3 samples, and the NaN it gives to the
    % *_dec features below it), the 1/e autocorrelation time used to coarse-grain in
    % 'rescale_tau' (EN_MSE), and the noise-gated AMI minimum of
    % BF_GetTau(y,'mi'). All series are made with BF_Random, so they are the
    % same in MATLAB and Python.

    methods (TestClassSetup)
        function runStartup(tc)
            try
                run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup.m'));
                pass = true;
            catch
                pass = false;
            end
            tc.fatalAssertTrue(pass, 'HCTSA failed to startup successfully.')
        end

        function clearSharedCaches(~)
            clear CO_AutoCorr CO_FirstMin BF_GetTau
        end
    end

    methods (Static)
        function y = noisyRandomWalk(seed)
            % a random walk observed with noise of 10 times its step size
            N = 1000;
            y = zscore(cumsum(0.1*BF_Random(N, seed, 'normal')) + BF_Random(N, seed + 100, 'normal'));
        end
    end

    methods (Test)
        %-------------------------------------------------------------------
        % 'decimate_ac1e': length floor
        %-------------------------------------------------------------------
        function decimateBelowFloorIsNaN(tc)
            % Fewer than 3 decimated samples gives NaN (a 2-sample series):
            y = zscore(BF_Random(2, 1, 'normal'));
            tc.verifyEqual(numel(1:BF_GetTau(y, 'ac1e'):numel(y)), 2);
            yDec = BF_PreProcess(y, 'decimate_ac1e');
            tc.verifyTrue(isscalar(yDec) && isnan(yDec));
            % ... and 3 samples are kept
            y = zscore(BF_Random(3, 1, 'normal'));
            tc.verifyEqual(BF_PreProcess(y, 'decimate_ac1e'), y(1:BF_GetTau(y, 'ac1e'):end));
        end

        function decimateShortResultIsKept(tc)
            % A random walk of 1000 samples decorrelates slowly and decimates to
            % fewer than 50 samples, which are kept (a 50-sample floor gave NaN):
            y = zscore(cumsum(BF_Random(1000, 1, 'normal')));
            tau = BF_GetTau(y, 'ac1e');
            tc.verifyLessThan(numel(1:tau:numel(y)), 50);
            tc.verifyEqual(BF_PreProcess(y, 'decimate_ac1e'), y(1:tau:end));
        end

        function decimateIsPlainSubsample(tc)
            % Prefixes of one AR(1) realization: every tau-th sample.
            u = filter(1, [1, -0.95], BF_Random(1500, 3, 'normal'));
            for i = 1:11
                y = zscore(u(1:400 + 100*i));
                tau = BF_GetTau(y, 'ac1e');
                tc.verifyEqual(BF_PreProcess(y, 'decimate_ac1e'), y(1:tau:end));
            end
        end

        function decimatedFeaturesAreNaNBelowFloor(tc)
            % Every registered master operation on the decimated series is NaN
            % (no error, and no finite value from the NaN input); a 2-sample
            % series is below the floor:
            x_z = zscore(BF_Random(2, 1, 'normal')); %#ok<NASGU> (used in eval)
            M = [];
            evalc('M = TS_ReadInputFile(''mops'', ''INP_mops_hctsa.txt'', false);');
            codes = M.Code(contains(M.Code, 'decimate_ac1e'));
            tc.verifyNotEmpty(codes);
            for i = 1:numel(codes)
                out = [];
                evalc('out = eval(codes{i});'); % (hide the operations' warnings)
                if isstruct(out)
                    vals = struct2cell(out);
                    tc.verifyTrue(all(cellfun(@(v) all(isnan(v)), vals)), codes{i});
                else
                    tc.verifyTrue(isscalar(out) && isnan(out), codes{i});
                end
            end
        end

        function guardedOperationsGiveNaNForNaNInput(tc)
            tc.verifyTrue(isnan(CO_TranslateShape(NaN, 'circle', 2.5, 'pts')));
            tc.verifyTrue(isnan(EX_MovingThreshold(NaN, 0.1, 0.1)));
            tc.verifyTrue(isnan(SB_BinaryStats(NaN, 'diff')));
            out = MD_Porta(NaN, 6);
            tc.verifyTrue(all(isnan(cell2mat(struct2cell(out)))));
        end

        %-------------------------------------------------------------------
        % 'rescale_tau': coarse-graining at the 1/e autocorrelation time
        %-------------------------------------------------------------------
        function rescaleTauUsesAC1eDelay(tc)
            y = zscore(filter(1, [1, -0.9], BF_Random(1000, 5, 'normal')));
            tau = BF_GetTau(y, 'ac1e');
            % (the first zero crossing of the ACF is much later, so the test
            % distinguishes the two)
            tc.verifyLessThan(tau, CO_FirstCrossing(y, 'ac', 0, 'discrete'));
            tc.verifyEqual(BF_PreProcess(y, 'rescale_tau'), mean(BF_MakeBuffer(y, tau), 2));
        end

        function multiscaleEntropyRescalesAtAC1e(tc)
            y = zscore(filter(1, [1, -0.9], BF_Random(1000, 5, 'normal')));
            tau = BF_GetTau(y, 'ac1e');
            yCG = zscore(mean(BF_MakeBuffer(y, tau), 2));
            outRescaled = EN_MSE(y, 1:10, 2, 0.15, 'rescale_tau');
            tc.verifyTrue(isstruct(outRescaled));
            tc.verifyEqual(outRescaled, EN_MSE(yCG, 1:10, 2, 0.15, ''));
        end

        %-------------------------------------------------------------------
        % BF_GetTau 'mi': noise-gated AMI minimum
        %-------------------------------------------------------------------
        function miIsNotOneForNoisyRandomWalk(tc)
            % The AMI curve of a noisy random walk is nearly flat at small lags,
            % so estimator noise can make the smoothed AMI rise from lag 1 to
            % lag 2 (seed 7, which gave tau = 1 without the noise gate). The
            % delay must not then collapse to 1 while the ACF decays slowly.
            for seed = 1:10
                y = SamplingRateGuardsTests.noisyRandomWalk(seed);
                tauAC = BF_GetTau(y, 'ac1e');
                if tauAC > 3
                    tc.verifyGreaterThan(BF_GetTau(y, 'mi'), 1, sprintf('seed %u', seed));
                end
            end
            y = SamplingRateGuardsTests.noisyRandomWalk(7);
            tc.verifyGreaterThan(BF_GetTau(y, 'ac1e'), 3);
            tc.verifyEqual(BF_GetTau(y, 'mi'), BF_GetTau(y, 'ac1e')); % (no significant minimum)
        end

        function miKeepsGenuineMinimumBelowAC1e(tc)
            % A slowly decorrelating AR(1) plus a period-20 oscillation: the
            % ACF only falls below 1/e after a few periods, but |ACF|, and so
            % the AMI, has a deep first minimum near half a period (lag 10).
            N = 1000;
            t = (1:N)';
            for seed = [2, 4, 8]
                slow = filter(1, [1, -0.98], BF_Random(N, seed, 'normal'));
                y = zscore(2*slow/std(slow) + sqrt(2)*sin(2*pi*t/20) + 0.3*BF_Random(N, seed + 50, 'normal'));
                tauAC = BF_GetTau(y, 'ac1e');
                tauMI = BF_GetTau(y, 'mi');
                tc.verifyGreaterThan(tauAC, 15, sprintf('seed %u', seed));
                tc.verifyLessThan(tauMI, tauAC, sprintf('seed %u', seed));
                tc.verifyGreaterThanOrEqual(tauMI, 6, sprintf('seed %u', seed));
                tc.verifyLessThanOrEqual(tauMI, 14, sprintf('seed %u', seed));
            end
        end

        function miOfOscillatorIsAtMostAC1e(tc)
            % A noisy AR(2) oscillation (period 40): its AMI minimum, at a
            % quarter period, is beyond the 1/e crossing, so 'mi' = 'ac1e'.
            T = 40;
            y = zscore(filter(1, [1, -2*0.97*cos(2*pi/T), 0.97^2], BF_Random(3000, 9, 'normal')));
            tauAC = BF_GetTau(y, 'ac1e');
            tc.verifyGreaterThan(tauAC, 3);
            tc.verifyEqual(BF_GetTau(y, 'mi'), tauAC);
        end
    end
end
