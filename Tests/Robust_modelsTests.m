classdef Robust_modelsTests < matlab.unittest.TestCase
    % Unit tests for model-fitting features that were made robust: bounded,
    % deterministic, and insensitive to tiny perturbations of the input.

    methods(TestClassSetup)
        function runStartup(testCase)
            try
                run("../startup.m")
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'HCTSA failed to startup successfully.')
        end
    end

    methods(Test)

        %-------------------------------------------------------------
        % CO_PartialAutoCorr
        %-------------------------------------------------------------
        function test_CO_PartialAutoCorr_AR1(testCase)
            % AR(1): partial autocorrelation at lag 1 is the coefficient, ~0 beyond
            rng(1);
            e = randn(5000,1); y = zeros(5000,1);
            for t = 2:5000, y(t) = 0.7*y(t-1) + e(t); end
            out = CO_PartialAutoCorr(y,5);
            testCase.verifyEqual(out.pac_1, 0.7, 'AbsTol', 0.03);
            testCase.verifyLessThan(abs(out.pac_2), 0.05);
            testCase.verifyLessThan(abs(out.pac_5), 0.05);
        end

        function test_CO_PartialAutoCorr_AgreesWithOLS(testCase)
            rng(2);
            y = filter(1,[1 -0.5 0.2],randn(1000,1));
            a = CO_PartialAutoCorr(y,10); b = CO_PartialAutoCorr(y,10,'ols');
            testCase.verifyEqual(cell2mat(struct2cell(a)), cell2mat(struct2cell(b)), 'AbsTol', 0.02);
        end

        function test_CO_PartialAutoCorr_BoundedAndStableOnSinusoid(testCase)
            % An exact sinusoid is perfectly predictable by an AR(2): the partial
            % autocorrelations beyond lag 2 are ~0 and do not depend on rounding noise
            t = (1:1000)';
            y = sin(2*pi*t/50);
            rng(3);
            a = cell2mat(struct2cell(CO_PartialAutoCorr(y,20)));
            b = cell2mat(struct2cell(CO_PartialAutoCorr(y + 1e-9*randn(1000,1),20)));
            testCase.verifyLessThanOrEqual(max(abs(a)), 1);
            testCase.verifyLessThan(max(abs(a(3:end))), 0.2);
            testCase.verifyEqual(a, b, 'AbsTol', 1e-3);
        end

        function test_CO_PartialAutoCorr_Constant(testCase)
            out = CO_PartialAutoCorr(ones(100,1),3);
            testCase.verifyEqual(cell2mat(struct2cell(out)), zeros(3,1));
        end

        %-------------------------------------------------------------
        % MF_hmm_Fit, MF_hmm_CompareNStates
        %-------------------------------------------------------------
        function test_MF_hmm_Fit_Deterministic(testCase)
            rng(4);
            y = zscore([randn(400,1); 3 + 0.5*randn(400,1)]);
            rng(1); a = MF_hmm_Fit(y,0.7,2);
            rng(99); b = MF_hmm_Fit(y,0.7,2);
            testCase.verifyEqual(a, b);
            testCase.verifyFalse(isfield(a,'nit'));
            testCase.verifyLessThan(a.Mu_1, a.Mu_2); % sorted means, two clear states
            testCase.verifyGreaterThan(a.Pmeandiag, 0.9);
        end

        function test_MF_hmm_Fit_FiniteOnQuantizedData(testCase)
            % Few distinct values: the variance floor keeps likelihoods finite
            rng(5);
            y = round(randn(500,1));
            out = MF_hmm_Fit(y,0.6,3);
            testCase.verifyTrue(all(structfun(@isfinite, out)));
            testCase.verifyGreaterThan(out.Cov, 0.01*var(y(1:300)) - 1e-12);
            out2 = MF_hmm_CompareNStates(y,0.6,2:4);
            testCase.verifyTrue(all(structfun(@isfinite, out2)));
        end

        function test_ZG_hmm_cl_NoUnderflow(testCase)
            % A test point far from every state gives a finite, large negative log likelihood
            [Mu,Cov,P,Pi] = deal([0;1], 0.01, [0.9 0.1;0.1 0.9], [0.5 0.5]);
            lik = ZG_hmm_cl([0; 50; 0], 3, 2, Mu, Cov, P, Pi);
            testCase.verifyTrue(isfinite(lik));
            testCase.verifyLessThan(lik, -1e4);
        end

        %-------------------------------------------------------------
        % MF_GP_*: noise floor
        %-------------------------------------------------------------
        function test_MF_GP_Hyperparameters_NoiseFloor(testCase)
            % A smooth noiseless series would push the noise to ~0 without a floor
            t = (1:200)';
            y = zscore(sin(2*pi*t/60) + 0.3*sin(2*pi*t/23));
            out = MF_GP_Hyperparameters(y,{'covSum',{'covSEiso','covNoise'}},1,200,'first');
            testCase.verifyGreaterThanOrEqual(out.logh3, log(0.01*std(y)) - 1e-8);
            testCase.verifyGreaterThanOrEqual(out.minS, 0.01*std(y) - 1e-8);
        end

        function test_MF_GP_FitAcross_ParameterizedCovariance(testCase)
            rng(6);
            y = zscore(cumsum(randn(200,1)));
            out = MF_GP_FitAcross(y,{'covSum',{{'covMaterniso',3},'covNoise'}},20);
            testCase.verifyTrue(isstruct(out));
            testCase.verifyFalse(isfield(out,'h_lonN'));
        end

        function test_MF_GP_LocalPrediction_SkipsConstantWindows(testCase)
            % Windows whose training data are constant cannot be standardized: they are
            % left out, so the ratios to the error bar stay moderate
            rng(10);
            y = randn(300,1); y(40:120) = 0.5;
            out = MF_GP_LocalPrediction(y,{'covSum',{'covSEiso','covNoise'}},10,3,20,'frombefore');
            testCase.verifyTrue(all(structfun(@isfinite, out)));
            testCase.verifyLessThan(out.maxabs_std, 1e4);
            testCase.verifyTrue(isnan(MF_GP_LocalPrediction(ones(300,1),{'covSum',{'covSEiso','covNoise'}},10,3,20,'frombefore')));
        end

        %-------------------------------------------------------------
        % WL_cwt
        %-------------------------------------------------------------
        function test_WL_cwt_ExactZerosGivesNaNGamma(testCase)
            y = [zeros(150,1); ones(150,1)];
            out = WL_cwt(y,'db3',8);
            testCase.verifyTrue(isnan(out.gam1));
            testCase.verifyTrue(isfinite(out.SC_h));
            rng(11);
            out2 = WL_cwt(randn(300,1),'db3',8);
            testCase.verifyTrue(isfinite(out2.gam1) && out2.gam1 > 0);
        end

        %-------------------------------------------------------------
        % PH_ForcePotential
        %-------------------------------------------------------------
        function test_PH_ForcePotential_MeanAbs(testCase)
            rng(7);
            y = zscore(randn(500,1));
            out = PH_ForcePotential(y,'dblwell',[1,0.5,0.2]);
            testCase.verifyFalse(isfield(out,'finaldev'));
            testCase.verifyGreaterThan(out.meanabs, 0);
            % mean |x| is at least |mean| and at most the largest excursion
            testCase.verifyGreaterThanOrEqual(out.meanabs, abs(out.mean) - 1e-12);
            testCase.verifyLessThanOrEqual(out.meanabs, out.range);
            out2 = PH_ForcePotential(y + 1e-12*randn(500,1),'dblwell',[1,0.5,0.2]);
            testCase.verifyEqual(out.meanabs, out2.meanabs, 'AbsTol', 1e-6);
        end

        %-------------------------------------------------------------
        % MF_arfit, MF_AR_arcov
        %-------------------------------------------------------------
        function test_MF_arfit_ExactSinusoidIsNaN(testCase)
            t = (1:500)';
            y = zscore(sin(2*pi*t/40));
            testCase.verifyTrue(isnan(MF_arfit(y,1,8,'sbc')));
            testCase.verifyTrue(isnan(MF_AR_arcov(y,4)));
        end

        function test_MF_arfit_AR2_Recovered(testCase)
            rng(8);
            y = zscore(filter(1,[1 -0.6 0.3],randn(2000,1)));
            o = MF_arfit(y,1,6,'sbc');
            testCase.verifyEqual(o.A1, 0.6, 'AbsTol', 0.06);
            testCase.verifyEqual(o.A2, -0.3, 'AbsTol', 0.06);
            a = MF_AR_arcov(y,2);
            testCase.verifyEqual(a.a2, -0.6, 'AbsTol', 0.06);
            testCase.verifyEqual(a.a3, 0.3, 'AbsTol', 0.06);
        end

        %-------------------------------------------------------------
        % DN_TailIndex: closed-form generalized Pareto shape
        %-------------------------------------------------------------
        function test_DN_TailIndex_GPDShape(testCase)
            rng(9);
            n = 20000;
            u = rand(n,1);
            yPareto = (1 - u).^(-0.5);   % tail index xi = 0.5 (above the median)
            yUnif = rand(n,1);            % xi = -1
            o1 = DN_TailIndex(yPareto, 0.05);
            o2 = DN_TailIndex(yUnif, 0.05);
            testCase.verifyEqual(o1.gpdUpper, 0.5, 'AbsTol', 0.12);
            testCase.verifyEqual(o2.gpdUpper, -1, 'AbsTol', 0.12);
            o1b = DN_TailIndex(yPareto, 0.05);
            testCase.verifyEqual(o1.gpdUpper, o1b.gpdUpper); % deterministic
        end

    end
end
