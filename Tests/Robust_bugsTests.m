classdef Robust_bugsTests < matlab.unittest.TestCase
    % Regression tests for bug fixes in DVV_surrogate, NL_PersistentHomology
    % and MF_CompareAR.

    methods (TestClassSetup)
        function setPaths(~)
            % Make the compiled ripser binary findable (as startup.m does)
            here = fileparts(fileparts(mfilename('fullpath')));
            setenv('PATH', [fullfile(here, 'Toolboxes', 'ripser', 'bin'), ':', getenv('PATH')]);
        end
    end

    methods (Test)
        function dvvSurrogateMatchesSpectrum(tc)
            % The surrogate must be rank-ordered by value: its amplitude spectrum then
            % approximates the original's (ranking by |s| leaves a ~100% mismatch for a
            % random walk)
            rng(5);
            x = zscore(cumsum(randn(400, 1)));
            Xs = DVV_surrogate(x, 2);
            A = abs(fft(x));
            for k = 1:2
                tc.verifyEqual(sort(Xs(:,k)), sort(x), 'AbsTol', 1e-12);
                tc.verifyLessThan(mean(abs(A - abs(fft(Xs(:,k))))) / mean(A), 0.1);
            end
        end

        function persistentHomologyH0IndependentOfMaxDim(tc)
            % totalPersistenceH0 must not pick up the dim-1 intervals when maxDim >= 1
            rng(2);
            t = (1:300)';
            y = zscore(sin(2*pi*t/25) + 0.1*randn(300, 1));
            o0 = NL_PersistentHomology(y, 5, 3, 0, 150);
            o1 = NL_PersistentHomology(y, 5, 3, 1, 150);
            tc.verifyEqual(o1.totalPersistenceH0, o0.totalPersistenceH0, 'AbsTol', 1e-10);
        end

        function compareARNaNOnShortSeries(tc)
            rng(4);
            o = MF_CompareAR(randn(20, 1), 1:10, 'all'); % highest order interpolates the data
            tc.verifyTrue(isnan(o));
            o = MF_CompareAR(randn(40, 1), 1:10, 0.5);   % training segment too short for these orders
            tc.verifyTrue(isnan(o));
            o = MF_CompareAR(randn(200, 1), 1:10, 0.5);
            tc.verifyGreaterThan(o.minv, 1e-3);
        end
    end
end
