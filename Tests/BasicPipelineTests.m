classdef BasicPipelineTests < matlab.unittest.TestCase
    
    methods(TestClassSetup)

        function runStartup(testCase)
            % run the startup script before executing the tests to ensure
            % that all paths and toolboxes are available
            try
                run("../startup.m")
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'HCTSA failed to startup successfully.')
        end

        function makeCopyofMat(testCase)
            % make a copy of the base .mat file
            duplicatedFileName = 'HCTSA_Bonn_EEG.mat';
            originalFileName = 'HCTSA_Bonn_EEG_OG.mat';

            % check original exists
            doesOGExist = exist(originalFileName, 'file') == 2;
            testCase.fatalAssertTrue(doesOGExist, 'Original HCTSA mat file not found.')
        
            % make a copy 
            copyfile(originalFileName, duplicatedFileName);
            % check for duplicated file
            doesDuplicateExist = exist(duplicatedFileName, 'file') == 2;
            testCase.fatalAssertTrue(doesDuplicateExist, 'Duplicated HCTSA mat file not found.')
        end          

        function checkRequireFiles(testCase)
            % check that all the files required for unit tests exist
            requiredFiles = {'INP_unit_test.mat', 'HCTSA_Bonn_EEG.mat', ...
                'HCTSA_catch22_expected.mat'};
            for i = 1:numel(requiredFiles)
                testCase.fatalAssertTrue(exist(requiredFiles{i}, 'file') == 2, ...
                    sprintf('Required file "%s" does not exist.', requiredFiles{i}));
            end
        end
      
    end
    
    methods(TestMethodSetup)
    end

    methods(TestClassTeardown)
        % clean up any remaining artifcats after unit tests have finished. 
         function cleanupTestClass(testCase)
                generatedFiles = {'HCTSA_catch22.mat', 'HCTSA_locDepFiltered.mat', ...
                    'HCTSA.mat', 'HCTSA_EEG_N.mat', 'HCTSA_Bonn_EEG.mat', ..., 
                    'HCTSA_Bonn_EEG_N.mat'};
                for i = 1:numel(generatedFiles)
                    if exist(generatedFiles{i}, 'file') == 2
                    delete(generatedFiles{i});
                    end
                end
    
                % close any opened figures
                close all;
           end
    end
    
    methods(Test)
        % Test methods
        
        function test_TS_Init(testCase)
            % before creating a new HCTSA.mat file, delete if already
            % exists to avoid user input prompt.
            expectedMatFile = 'HCTSA.mat';
            if exist(expectedMatFile, 'file')
                delete(expectedMatFile)
            end

            % specify the input .mat file
            inputFile = 'INP_unit_test.mat';

            % call the TS_Init function in non-interactive mode
            try
                TS_Init(inputFile, 'hctsa', false);
                pass = true;
            catch
                % catch any exceptions raised
                pass = false;
            end
             % check whether TS_Compute executes successfully
            testCase.fatalAssertTrue(pass, 'TS_Init did not execute sucessfully.')
                

            % check that the HCTSA.mat file was created
            fileExists = exist(expectedMatFile, 'file') == 2;
            testCase.fatalAssertTrue(fileExists, 'HCTSA.mat file was not created');
           
        end

        function test_TS_CalculateFeatureVector(testCase)
            % test that the function executes successfully
            x = rand(500, 1); % random time series
            try
                featVector = TS_CalculateFeatureVector(x,false);
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_CalculateFeatureVector did not execute sucessfully.')
            % inspect output
            num_ts = size(featVector, 2);
            testCase.verifyEqual(num_ts, 1, 'Feature vector output shape unexpected (not a vector).')
        end

        function test_TS_ComputeCatch22(testCase)
            % test catch22 subset on 31 sample time series and compare to the expected output.
            % check if HCTSA_catch22.mat exists and delete if true. 
            % assumes that the TS_Init function works as expected from
            % previous unit test. 

            if exist('HCTSA_catch22.mat', 'file')
                delete('HCTSA_catch22.mat')
            end
            TS_Init('INP_unit_test.mat','catch22',false,'HCTSA_catch22.mat');
            try
                TS_Compute(false, [], [], 'missing', 'HCTSA_catch22.mat');
                pass = true;
            catch
                % catch any exceptions raised
                pass = false;
            end
             % check whether TS_Compute executes successfully
            testCase.fatalAssertTrue(pass, 'TS_Compute did not execute sucessfully for catch22 set.')

            % compare to expected output
            actual_output = load('HCTSA_catch22.mat', 'TS_DataMat').TS_DataMat;
            expected_output = load('HCTSA_catch22_expected.mat', 'TS_DataMat').TS_DataMat;
            testCase.verifyEqual(actual_output, expected_output, 'Expected output != actual output, catch22.', "RelTol",0.1)

        end

        function test_TS_ComputeHCTSA(testCase)
            % test the entire feature set on 1 sample time series
            try
                TS_Compute(false, [], [], 'missing', '', false );
                pass = true;
            catch
                % catch any exceptions raised
                pass = false;
            end
            % check whether TS_Compute executes successfully
            testCase.fatalAssertTrue(pass, 'TS_Compute did not execute sucessfully for hctsa set.')

            % run basic checks on the output
            actual_output = load("HCTSA.mat", "TS_DataMat").TS_DataMat;

            % (1) Check that the output contains numerics (not NaNs) -
            % which would indicate that features values weren't computed or
            % filled in.
            testCase.verifyFalse(all(isnan(actual_output), 'all'), 'TS_DataMat contains all NaNs. Features not computed.')
            

            % (2) Check the size of the output
            
            num_ts = size(actual_output, 1);
            testCase.verifyEqual(num_ts, 1, 'Expected 1 time series.')


            % (3) Compare actual output to expected output.
            % load in benchmark values
            %expected_output = load('HCTSA_expected.mat', 'TS_DataMat').TS_DataMat
            %testCase.verifyEqual(actual_output, expected_output, 'Expected output not equal to actual output on test TS.', "RelTol",0.1)
            
        end

        function test_TS_InspectQuality(testCase)
            % simple check of whether or not the function runs. 
            try
                out = TS_InspectQuality('summary', 'HCTSA_Bonn_EEG.mat');
                close()
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_InspectQuality did not execute sucessfully.')
            % generate summary report of errors

        end

        function test_TS_Normalize(testCase)
            % check whether TS_Normalize executes successfully (w. default)
            expectedOutputMatFile = 'HCTSA_Bonn_EEG_N.mat';
            if exist(expectedOutputMatFile, 'file')
                delete(expectedOutputMatFile);
            end
            try
                TS_Normalize('', [], 'HCTSA_Bonn_EEG.mat');
                pass = true;
            catch
                 pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_Compute did not execute sucessfully for hctsa set.')

            % check for generation of the normalized .mat file
            fileExists = exist(expectedOutputMatFile, 'file') == 2;
            testCase.fatalAssertTrue(fileExists, 'Normalized .mat file was not created');

        end

        function test_TS_PlotTimeSeries(testCase)
            % basic check - does the function run
            try
                TS_PlotTimeSeries('HCTSA_Bonn_EEG.mat')
                close()
                pass = true;
            catch
                pass = false;
            end

            testCase.fatalAssertTrue(pass, 'TS_PlotTimeSeries did not execute sucessfully.')

        end

        function test_TS_PlotDataMatrix(testCase)
            % basic check - does the function run
            try
                TS_PlotDataMatrix('whatData','HCTSA_Bonn_EEG_N.mat') % uses norm. data
                close()
                pass = true;
            catch
                pass = false;
            end

            testCase.fatalAssertTrue(pass, 'TS_PlotDataMatrix did not execute sucessfully.')

        end

        function test_TS_Subset(testCase)
            %  basic check - does the function run without errors? 
            newFilteredMatName = 'HCTSA_locDepFiltered.mat';
            % delte new filtered mat if already exists before unit test
            if exist(newFilteredMatName, 'file')
                delete(newFilteredMatName)
            end
            try
                % filter out location dependent features
                [IDs_locDep,IDs_notLocDep] = TS_GetIDs('locdep','HCTSA_Bonn_EEG.mat','ops','Keywords');
                TS_Subset('HCTSA_Bonn_EEG.mat',[],IDs_notLocDep,true,newFilteredMatName);
                pass = true;
            catch
                % 
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_Subset did not execute sucessfully.');

            % check that the new filtered .mat file exists
            doesFileExist = exist(newFilteredMatName, 'file') == 2;
            testCase.verifyTrue(doesFileExist, 'Filtered .mat file was not created.')

            % delete artifacts
            delete(newFilteredMatName)
        end

        function test_TS_PlotLowDim(testCase)
            % basic check - does the function run without errors
            try
                TS_PlotLowDim('HCTSA_Bonn_EEG_N.mat', 'pca')
                close()
                TS_PlotLowDim('HCTSA_Bonn_EEG_N.mat', 'tsne')
                close()
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_PlotLowDim did not execute sucessfully.');
        end

        function test_TS_LabelGroups(testCase)
            % basic check - does the function run without errors
            try
                TS_LabelGroups('HCTSA_Bonn_EEG_N.mat', {'eyesOpen', 'seizure'});
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_LabelGroups did not execute sucessfully.');
        end

        function test_GiveMeDefaultClassificationParams(testCase)
            % test the default classification params when using different
            % labels
            % (2 class problem)
            TS_LabelGroups('HCTSA_Bonn_EEG_N.mat', {'eyesOpen', 'eyesClosed'});
            try
                cfnParams = GiveMeDefaultClassificationParams('HCTSA_Bonn_EEG_N.mat');
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'GiveMeDefaultClassificationParams did not execute sucessfully.');
            % check fields of struct
            isTwoClasses = cfnParams.numClasses == 2;
            testCase.verifyTrue(isTwoClasses, 'Expected 2 class problem for 2 labels.');

            % (5 class problem)
            TS_LabelGroups('HCTSA_Bonn_EEG_N.mat', {'eyesOpen', 'seizure', 'eyesClosed', 'epileptogenic', 'hippocampus'});
            cfnParams = GiveMeDefaultClassificationParams('HCTSA_Bonn_EEG_N.mat');
            isFiveClasses = cfnParams.numClasses == 5;
            testCase.verifyTrue(isFiveClasses, 'Expected 5 class problem for 5 labels.');

            % check setting of fields
            cfnParams = GiveMeDefaultClassificationParams('HCTSA_Bonn_EEG_N.mat');
            cfnParams.numRepeats = 10;
            isNumRepeatsCorrect = cfnParams.numRepeats == 10;
            testCase.verifyTrue(isNumRepeatsCorrect, 'cfnParams num repeats not correctly set.');
            cfnParams.numFolds = 5;
            isNumFoldsCorrect = cfnParams.numFolds == 5;
            testCase.verifyTrue(isNumFoldsCorrect, 'cfnParams num folds not correctly set.');
        end

        function test_TS_Classify(testCase)
            % test default settings
            TS_LabelGroups('HCTSA_Bonn_EEG_N.mat', {'eyesOpen', 'seizure', 'eyesClosed', 'epileptogenic', 'hippocampus'});
            try
                TS_Classify('HCTSA_Bonn_EEG_N.mat', struct(), '', doPlot=false);
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_Classify did not execute sucessfully with default params.');

            % test setting of num nulls
            TS_LabelGroups('HCTSA_Bonn_EEG_N.mat', {'epileptogenic','hippocampus'});
            numNulls = 10;
            try
                TS_Classify('HCTSA_Bonn_EEG_N.mat', struct(), numNulls, doPlot=false);
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_Classify did not execute sucessfully with num nulls setting.');
        end

        function test_TS_CompareFeatureSets(testCase)
            % basic checks - does the function run without error
            TS_LabelGroups('HCTSA_Bonn_EEG_N.mat', {'epileptogenic','hippocampus'});
            cfnParams = GiveMeDefaultClassificationParams('HCTSA_Bonn_EEG_N.mat');
            % use params for fast evaluation/basic checks for successful
            % execution
            cfnParams.numFolds = 2;
            cfnParams.numRepeats = 1;
            try
                TS_CompareFeatureSets('HCTSA_Bonn_EEG_N.mat',cfnParams);
                close();
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_CompareFeatureSets did not execute sucessfully.');
        end

        function test_TS_ClassifyLowDim(testCase)
            % basic checks - does the function run without error
            TS_LabelGroups('HCTSA_Bonn_EEG_N.mat', {'epileptogenic','hippocampus'});
            cfnParams = GiveMeDefaultClassificationParams('HCTSA_Bonn_EEG_N.mat');
            cfnParams.numFolds = 2;
            cfnParams.numRepeats = 1;
            try
                TS_ClassifyLowDim('HCTSA_Bonn_EEG_N.mat', cfnParams, 5, false);
                close();
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_ClassifyLowDim did not execute sucessfully.');
        end

        function test_TS_TopFeatures(testCase)
            % basic checks - does the function run without error
            TS_LabelGroups('HCTSA_Bonn_EEG_N.mat', {'epileptogenic','hippocampus'})
            try
                TS_TopFeatures('HCTSA_Bonn_EEG_N.mat','classification', struct(), 'whatPlots',{'histogram'});
                close();
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_TopFeatures did not execute sucessfully.');
        end

        function test_TS_FeatureSummary(testCase)
            % basic checks - does the function run without error
            TS_LabelGroups('HCTSA_Bonn_EEG_N.mat', {'epileptogenic','hippocampus'})
            try
                TS_FeatureSummary(95, 'HCTSA_Bonn_EEG_N.mat', true, false);
                close();
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_FeatureSummary did not execute sucessfully.');
        end

        function test_TS_SimSearch(testCase)
            % basic checks - does the function run without error
            try
                TS_SimSearch('whatData', 'HCTSA_Bonn_EEG_N.mat');
                close();
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_SimSearch did not execute sucessfully.');

            % test network viz of op 30 neighbours
            try
                % generate network visualisation for operation 30
                TS_SimSearch('whatData', 'HCTSA_Bonn_EEG_N.mat', 'targetID', 30,'whatPlots',{'network'}, 'tsOrOps', 'ops');
                close();
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_SimSearch network visualisation did not execute sucessfully.');
        end

        function test_TS_SingleFeature(testCase)
            % basic checks - does the function run without error
            opID = 500; % choose a random operation
            makeViolin = false;
            try
                TS_SingleFeature('HCTSA_Bonn_EEG_N.mat', opID, makeViolin);
                close(); % close out the opened figure
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_SingleFeature did not execute sucessfully.');
        end

    end
    
end
