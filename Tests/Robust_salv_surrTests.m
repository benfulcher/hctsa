classdef Robust_salv_surrTests < matlab.unittest.TestCase
    % SD_Surrogates: the meannumsurr output (mean of the statistic's numerator over the
    % surrogates).

    methods (Test)
        function meannumsurrIsMeanOfNumerator(testCase)
            % meannumsurr is the mean of CO_TC3/CO_trev's numerator over the same surrogates
            y = zscore(BF_Random(300, 7, 'normal'));
            y = zscore(y.^2 + 0.5 * [0; y(1:end-1)]);
            for surrMethod = 1:3
                for fn = {'tc3', 'trev'}
                    out = SD_Surrogates(y, 2, 20, surrMethod, fn{1}, 3);
                    names = {'RP', 'AAFT', 'RandPerm'};
                    surr = SD_MakeSurrogates(y, names{surrMethod}, 20, [], 3);
                    num = zeros(20, 1);
                    for i = 1:20
                        if strcmp(fn{1}, 'tc3')
                            num(i) = CO_TC3(surr(:, i), 2).num;
                        else
                            num(i) = CO_trev(surr(:, i), 2).num;
                        end
                    end
                    testCase.verifyEqual(out.meannumsurr, mean(num), 'AbsTol', 1e-12);
                end
            end
        end

        function aaftNumeratorTracksSkewedLinearProcess(testCase)
            % A static skewing transform of a linear Gaussian process gives a positive
            % third-order moment in its AAFT surrogates; a Gaussian process gives ~0
            e = BF_Random(2000, 11, 'normal');
            g = filter(1, [1, -0.8], e);
            gauss = zscore(g);
            skewed = zscore(exp(g / std(g)));
            oG = SD_Surrogates(gauss, 1, 100, 2, 'tc3', 0);
            oS = SD_Surrogates(skewed, 1, 100, 2, 'tc3', 0);
            testCase.verifyGreaterThan(oS.meannumsurr, 0.5);
            testCase.verifyLessThan(abs(oG.meannumsurr), 0.1);
            % reproducible with the same seed
            oS2 = SD_Surrogates(skewed, 1, 100, 2, 'tc3', 0);
            testCase.verifyEqual(oS2.meannumsurr, oS.meannumsurr);
        end
    end
end
