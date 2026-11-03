classdef AnalysisFunctionsTests < matlab.unittest.TestCase
    % Numerical behavior of the analysis helper functions that do not need computed
    % features: BF_NaNCov on features with very large values, and the ordering of
    % matches in TS_SimSearch.

    methods (TestClassSetup)
        function setUp(testCase)
            run(fullfile(fileparts(fileparts(mfilename('fullpath'))),'startup.m'));
            oldState = get(groot,'defaultFigureVisible');
            set(groot,'defaultFigureVisible','off');
            testCase.addTeardown(@() set(groot,'defaultFigureVisible',oldState));
        end
    end

    methods (Test)
        function nanCovHandlesValuesBeyondSinglePrecision(testCase)
            % Values above ~3.4e38 are not representable in single precision, nor are
            % their products, which used to make every correlation with them infinite
            rng(1);
            X = randn(200,8);
            X(:,2) = 1e100*X(:,2) + 3e100*X(:,1);
            X(:,3) = 1e200*X(:,3); % (its sum of squares is beyond double precision)
            X(7,[1,2,5]) = NaN;
            R = double(BF_NaNCov(X,1,0));
            testCase.verifyTrue(all(isfinite(R(:))));
            testCase.verifyEqual(diag(R),ones(8,1),'AbsTol',1e-3);
            refR = corr(X./max(abs(X)),'Rows','pairwise'); % (corr itself overflows on the huge columns)
            testCase.verifyEqual(R,refR,'AbsTol',2e-2);
            % ... also for the distance form used for 'corr_fast':
            D = BF_NaNCov(X,1,1);
            testCase.verifyTrue(all(isfinite(D(:))));
        end

        function nanCovHandlesValuesBeyondDoublePrecisionSquares(testCase)
            % With and without NaNs, for features around 1e200 (variances overflow double
            % precision unless the columns are scaled first)
            rng(5);
            X = randn(150,5);
            X(:,2) = 1e200*X(:,2);
            X(:,4) = 3e250*X(:,4) + 1e250*X(:,1);
            for hasNaN = [false,true]
                Y = X;
                if hasNaN, Y(9,[1,2]) = NaN; end
                R = double(BF_NaNCov(Y,1,0));
                testCase.verifyTrue(all(isfinite(R(:))),sprintf('NaNs: %d',hasNaN));
                testCase.verifyEqual(R,corr(Y./max(abs(Y)),'Rows','pairwise'),'AbsTol',2e-2);
                % Covariances are finite where they are representable:
                C = double(BF_NaNCov(Y,0,0));
                rows = ~isnan(Y(:,2));
                refC = cov(Y(rows,2)/1e200,Y(rows,5));
                testCase.verifyEqual(C(2,5),refC(1,2)*1e200,'RelTol',5e-2);
            end
        end

        function nanCovIsScaleInvariantForCorrelations(testCase)
            rng(2);
            X = randn(100,6);
            X(3,2) = NaN;
            Xbig = X.*[1,1e50,1,1e-50,1e150,1];
            testCase.verifyEqual(double(BF_NaNCov(Xbig,1,0)),double(BF_NaNCov(X,1,0)),'AbsTol',1e-5);
        end

        function nanCovMatchesCovForWellScaledData(testCase)
            rng(3);
            X = randn(300,20);
            X(5,3) = NaN;
            C = BF_NaNCov(X,0,0);
            refC = cov(X,'omitrows'); % (approximately the same; the means differ slightly)
            testCase.verifyEqual(C,refC,'AbsTol',0.05);
            % Without NaNs it is exactly cov:
            testCase.verifyEqual(BF_NaNCov(X(:,[1:2,4:20]),0,0),cov(X(:,[1:2,4:20])));
        end

        function simSearchPutsTargetFirstAndNotOverflowedFeature(testCase)
            % On raw data with a feature above 3.4e38, 'corr_fast' used to give that
            % feature distance -Inf to every target, which tied with the target's own
            % sort key and replaced it as the first match
            tmp = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(tmp.Folder));
            here = fileparts(mfilename('fullpath'));
            matFile = fullfile(here,'HCTSA_Bonn_EEG_OG.mat');
            if ~isfile(matFile)
                testCase.assumeFail('HCTSA_Bonn_EEG_OG.mat not found');
            end
            S = load(matFile,'Operations','TimeSeries');
            % A subset of time series and operations, with all values good:
            hasData = load(matFile,'TS_DataMat','TS_Quality');
            goodOps = find(all(hasData.TS_Quality(1:40,:)==0,1));
            opInd = goodOps(1:60);
            tsIDs = S.TimeSeries.ID(1:40);
            opIDs = S.Operations.ID(opInd);
            TS_Subset(matFile,tsIDs,opIDs,true,'HCTSA.mat');

            load('HCTSA.mat','Operations');
            % Replace the values by independent random features, except for a near-copy
            % of the target, one huge feature, and one with a missing value (so that
            % the NaN-tolerant path is used):
            rng(4);
            TS_DataMat = randn(40,60);
            hugeCol = 2; targetCol = 5; copyCol = 6; nanCol = 8; % (huge before the target: a stable sort lists ties in this order)
            TS_DataMat(:,copyCol) = TS_DataMat(:,targetCol) + 0.1*randn(40,1);
            TS_DataMat(1,targetCol) = 5; TS_DataMat(1,copyCol) = 5; % (positive, to match the sign below)
            TS_DataMat(1,hugeCol) = 1e40; % (an outlier, so only this row overflows in single precision)
            TS_DataMat(3,nanCol) = NaN;
            save('HCTSA.mat','TS_DataMat','-append');

            output = evalc(['TS_SimSearch(''targetID'',Operations.ID(targetCol),''tsOrOps'',''ops'',' ...
                            '''whatData'',''HCTSA.mat'',''whatDistMetric'',''corr_fast'',''whatPlots'',{},' ...
                            '''numNeighbors'',3)']);
            tokens = regexp(output,'^1\. \[(\d+)\]','tokens','lineanchors','once');
            testCase.verifyNotEmpty(tokens,'No first match was listed');
            testCase.verifyEqual(str2double(tokens{1}),Operations.ID(copyCol));
            % The target itself is never listed as a match:
            testCase.verifyEmpty(regexp(output,sprintf('\\d\\. \\[%u\\]',Operations.ID(targetCol)),'once'));
        end
    end
end
