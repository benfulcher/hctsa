classdef BasicPipelineTests < matlab.unittest.TestCase

    properties
        % log outputs
        logFileName = 'summary_log.txt';
        failedOpFileName = 'failed_operations_benchmark.txt';
        footer = [repmat('=', 1, 62), '\n'];
        % benchmarking files
        precomputedMatNameOG = 'HCTSA_Bonn_EEG_OG.mat';
        precomputedMatName = 'HCTSA_Bonn_EEG.mat';
        precomputedMatNormName = 'HCTSA_Bonn_EEG_N.mat'
        precomputedCatch22Mat = 'HCTSA_catch22_expected.mat';
        timeSeriesINP = 'INP_unit_test.mat';
        requiredFiles;
    end

    methods(TestClassSetup)

        function runStartup(testCase)
            try
                run("../startup.m")
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'HCTSA failed to startup successfully.')
            % set the required files for benchmarking
            testCase.requiredFiles = {testCase.timeSeriesINP, testCase.precomputedMatNameOG, testCase.precomputedCatch22Mat};
            
        end

        function makeCopyofMat(testCase)
            % make a copy of the base .mat file
            duplicatedFileName = testCase.precomputedMatName;
            originalFileName = testCase.precomputedMatNameOG;

            doesOGExist = exist(originalFileName, 'file') == 2;
            testCase.fatalAssertTrue(doesOGExist, 'Original HCTSA mat file not found.')
        
            copyfile(originalFileName, duplicatedFileName);
            doesDuplicateExist = exist(duplicatedFileName, 'file') == 2;
            testCase.fatalAssertTrue(doesDuplicateExist, 'Duplicated HCTSA .mat file not found.')
        end          

        function checkRequiredFiles(testCase)
            % check that all the files required for unit tests exist
            req_files = testCase.requiredFiles;
            for i = 1:numel(req_files)
                testCase.fatalAssertTrue(exist(req_files{i}, 'file') == 2, ...
                    sprintf('Required file "%s" does not exist.', req_files{i}));
            end
        end

    end
    
    methods(TestClassTeardown)
        % clean up any remaining artifcats after unit tests have finished. 
         function cleanupTestClass(testCase)
                allMatFiles = dir('*.mat');
                allMatFiles = {allMatFiles.name};
                filesToDelete = setdiff(allMatFiles, testCase.requiredFiles);
                for i = 1:numel(filesToDelete)
                    if exist(filesToDelete{i}, 'file') == 2
                        delete(filesToDelete{i});
                    end
                end
                % tidy up
                close all;
                % print summary table
                testCase.printLogFile();
                % delete log files
                %delete(testCase.logFileName);
                %delete(testCase.failedOpFileName);
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
            inputFile = testCase.timeSeriesINP;

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
            % test three time series of varying lengths
            lengths = [500, 1000, 5000];
            f = 10;
            fs = 100;
            header = sprintf('\n  -------HCTSA Execution Time Benchmarks-------   \n');
            testCase.writeToLog(header);

            for i = 1:length(lengths)
                N = lengths(i);
                t = (0:N-1)/fs;
                x = 1.0 * sin(2*pi*f*t);
                x_noisy = x + 0.1 * randn(size(x));
                try
                    startTime = tic;
                    featVector = TS_CalculateFeatureVector(x_noisy');
                    executionTime = toc(startTime);
                    pass = true;
                catch
                    pass = false;
                end
                message = sprintf("Noisy sinusoid of length %d samples, time taken: %.4f seconds", lengths(i), executionTime);
                testCase.writeToLog(message);
                testCase.fatalAssertTrue(pass, 'TS_CalculateFeatureVector did not execute sucessfully.')
                num_ts = size(featVector, 2);
                testCase.verifyEqual(num_ts, 1, 'Feature vector output shape unexpected (not a vector).')
            end
        end

        function test_TS_ComputeCatch22(testCase)
            % test catch22 subset on sample time series and compare to the expected output.
            % check if HCTSA_catch22.mat exists and delete if true. 
            % assumes that the TS_Init function works as expected from
            % previous unit test. 
            matName = 'HCTSA_catch22.mat';

            if exist(matName, 'file')
                delete(matName)
            end
            TS_Init(testCase.timeSeriesINP,'catch22',false,matName);
            try
                TS_Compute(false, [], [], 'missing', matName);
                pass = true;
            catch
                % catch any exceptions raised
                pass = false;
            end
             % check whether TS_Compute executes successfully
            testCase.fatalAssertTrue(pass, 'TS_Compute did not execute sucessfully for catch22 set.')

            % compare to expected output
            actual_output = load(matName, 'TS_DataMat').TS_DataMat;
            expected_output = load(testCase.precomputedCatch22Mat, 'TS_DataMat').TS_DataMat;
            testCase.verifyEqual(actual_output, expected_output, 'Expected output != actual output, catch22.', "RelTol", 0.1)

        end

        function test_TS_ComputeHCTSA(testCase)
            % test the entire feature set on 1 sample time series
            try
                TS_Compute(false, [], [], 'missing', '', false);
                pass = true;
            catch
                % catch any exceptions raised
                pass = false;
            end
            % check whether TS_Compute executes successfully
            testCase.fatalAssertTrue(pass, 'TS_Compute did not execute sucessfully for hctsa set.')

            % run basic checks on the output
            actual_output = load("HCTSA.mat", "TS_DataMat").TS_DataMat;
            quality_labels = load("HCTSA.mat", "TS_Quality").TS_Quality;

            % (1) Check that the output contains numerics (not NaNs)
            testCase.verifyFalse(all(isnan(actual_output), 'all'), 'TS_DataMat contains all NaNs. Features not computed.')
            
            % (2) Check the quality labels
            quality_labels_counts = zeros(1, 8);
            for i = 0:7
                quality_labels_counts(i+1) = length(find(quality_labels == i));
            end


            header = sprintf('\n  -------Noisy Sinusoid Benchmark Feature Outputs-------   \n');
            testCase.writeToLog(header)
            quality_labels_names = {'successfully computed features','fatal errors','NaNs','positive infinities',...
                'negative infinities','imaginary values','empty outputs','non-existent features'};
            
            ops_all = load("HCTSA.mat", "Operations").Operations;
            num_ops_total = size(ops_all, 1);
            for i = 1:length(quality_labels_counts)
                percentage_op = (quality_labels_counts(i) / num_ops_total) * 100.0;
               message = sprintf("%s: %d [%.2f%%]", quality_labels_names{i}, quality_labels_counts(i), percentage_op);
                testCase.writeToLog(message);
            end
            testCase.writeToLog(testCase.footer)

            % (2) Check the size of the output
            num_ts = size(actual_output, 1);
            testCase.verifyEqual(num_ts, 1, 'Expected 1 time series.')

            % (3) Write all fatal error operations
            quality_labels = load("HCTSA.mat", "TS_Quality").TS_Quality;
            fatal_error_ops_ids = quality_labels == 1;
            fatal_error_ops_table = ops_all(fatal_error_ops_ids, :);
            writetable(fatal_error_ops_table, testCase.failedOpFileName);
        end

        function test_TS_InspectQuality(testCase)
            % simple check of whether or not the function runs. 
            try
                TS_InspectQuality('summary', testCase.precomputedMatName);
                close()
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_InspectQuality did not execute sucessfully.')

        end

        function test_TS_Normalize(testCase)
            % check whether TS_Normalize executes successfully (w. default)
            expectedNormMatFile = testCase.precomputedMatNormName;

            if exist(expectedNormMatFile, 'file')
                delete(expectedNormMatFile);
            end
            try
                TS_Normalize('', [], testCase.precomputedMatName);
                pass = true;
            catch
                 pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_Compute did not execute sucessfully for hctsa set.')

            % check for generation of the normalized .mat file
            fileExists = exist(expectedNormMatFile, 'file') == 2;
            testCase.fatalAssertTrue(fileExists, 'Normalized .mat file was not created');

        end

        function test_TS_PlotTimeSeries(testCase)
            % basic check - does the function run
            try
                TS_PlotTimeSeries(testCase.precomputedMatName)
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
                TS_PlotDataMatrix('whatData', testCase.precomputedMatNormName) % uses norm data
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
                [~,IDs_notLocDep] = TS_GetIDs('locdep',testCase.precomputedMatName,'ops','Keywords');
                TS_Subset(testCase.precomputedMatName,[],IDs_notLocDep,true,newFilteredMatName);
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
                TS_PlotLowDim(testCase.precomputedMatNormName, 'pca')
                close()
                TS_PlotLowDim(testCase.precomputedMatNormName, 'tsne')
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
                TS_LabelGroups(testCase.precomputedMatNormName, {'eyesOpen', 'seizure'});
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_LabelGroups did not execute sucessfully.');
        end

        function test_GiveMeDefaultClassificationParams(testCase)
            % test the default classification params when using different
            % labels (2 class problem)
            TS_LabelGroups(testCase.precomputedMatNormName, {'eyesOpen', 'eyesClosed'});
            try
                cfnParams = GiveMeDefaultClassificationParams(testCase.precomputedMatNormName);
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'GiveMeDefaultClassificationParams did not execute sucessfully.');
            % check fields of struct
            isTwoClasses = cfnParams.numClasses == 2;
            testCase.verifyTrue(isTwoClasses, 'Expected 2 class problem for 2 labels.');

            % (5 class problem)
            TS_LabelGroups(testCase.precomputedMatNormName, {'eyesOpen', 'seizure', 'eyesClosed', 'epileptogenic', 'hippocampus'});
            cfnParams = GiveMeDefaultClassificationParams(testCase.precomputedMatNormName);
            isFiveClasses = cfnParams.numClasses == 5;
            testCase.verifyTrue(isFiveClasses, 'Expected 5 class problem for 5 labels.');

            % check setting of fields
            cfnParams = GiveMeDefaultClassificationParams(testCase.precomputedMatNormName);
            cfnParams.numRepeats = 10;
            isNumRepeatsCorrect = cfnParams.numRepeats == 10;
            testCase.verifyTrue(isNumRepeatsCorrect, 'cfnParams num repeats not correctly set.');
            cfnParams.numFolds = 5;
            isNumFoldsCorrect = cfnParams.numFolds == 5;
            testCase.verifyTrue(isNumFoldsCorrect, 'cfnParams num folds not correctly set.');
        end

        function test_TS_Classify(testCase)
            % test default settings
            TS_LabelGroups(testCase.precomputedMatNormName, {'eyesOpen', 'seizure', 'eyesClosed', 'epileptogenic', 'hippocampus'});
            try
                TS_Classify(testCase.precomputedMatNormName, struct(), '', doPlot=false);
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_Classify did not execute sucessfully with default params.');

            % test setting of num nulls
            TS_LabelGroups(testCase.precomputedMatNormName, {'epileptogenic','hippocampus'});
            numNulls = 10;
            try
                TS_Classify(testCase.precomputedMatNormName, struct(), numNulls, doPlot=false);
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_Classify did not execute sucessfully with num nulls setting.');
        end

        function test_TS_CompareFeatureSets(testCase)
            % basic checks - does the function run without error
            TS_LabelGroups(testCase.precomputedMatNormName, {'epileptogenic','hippocampus'});
            cfnParams = GiveMeDefaultClassificationParams(testCase.precomputedMatNormName);
            % use params for fast evaluation/basic checks for successful
            % execution
            cfnParams.numFolds = 2;
            cfnParams.numRepeats = 1;
            try
                TS_CompareFeatureSets(testCase.precomputedMatNormName,cfnParams);
                close();
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_CompareFeatureSets did not execute sucessfully.');
        end

        function test_TS_ClassifyLowDim(testCase)
            % basic checks - does the function run without error
            TS_LabelGroups(testCase.precomputedMatNormName, {'epileptogenic','hippocampus'});
            cfnParams = GiveMeDefaultClassificationParams(testCase.precomputedMatNormName);
            cfnParams.numFolds = 2;
            cfnParams.numRepeats = 1;
            try
                TS_ClassifyLowDim(testCase.precomputedMatNormName, cfnParams, 5, false);
                close();
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_ClassifyLowDim did not execute sucessfully.');
        end

        function test_TS_TopFeatures(testCase)
            % basic checks - does the function run without error
            TS_LabelGroups(testCase.precomputedMatNormName, {'epileptogenic','hippocampus'})
            try
                TS_TopFeatures(testCase.precomputedMatNormName,'classification', struct(), 'whatPlots',{'histogram'});
                close();
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_TopFeatures did not execute sucessfully.');
        end

        function test_TS_FeatureSummary(testCase)
            % basic checks - does the function run without error
            TS_LabelGroups(testCase.precomputedMatNormName, {'epileptogenic','hippocampus'})
            try
                TS_FeatureSummary(95, testCase.precomputedMatNormName, true, false);
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
                TS_SimSearch('whatData', testCase.precomputedMatNormName);
                close();
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_SimSearch did not execute sucessfully.');

            % test network viz of op 30 neighbours
            try
                % generate network visualisation for operation 30
                TS_SimSearch('whatData', testCase.precomputedMatNormName, 'targetID', 30,'whatPlots',{'network'}, 'tsOrOps', 'ops');
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
                TS_SingleFeature(testCase.precomputedMatNormName, opID, makeViolin);
                close(); % close out the opened figure
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'TS_SingleFeature did not execute sucessfully.');
        end

    end

    methods(Access = private)
        function writeToLog(testCase, message)
            logFile = fopen(testCase.logFileName, 'a');
            fprintf(logFile, '%s\n', message);
            fclose(logFile);
        end

        function printLogFile(testCase)
            logFile = fopen(testCase.logFileName, 'r');
            logData = textscan(logFile, '%s', 'Delimiter', '\n');
            fclose(logFile);
            fprintf('\n==============================================================\n');
            fprintf('                    Unit Test Log                      \n');
            fprintf('==============================================================\n\n');
            % print the contents of the custom log file
            for i = 1:length(logData{1})
                fprintf('%s\n', logData{1}{i});
            end
            % print the table of fatal error operations
            fprintf('\n==============================================================\n');
            fprintf('      Noisy Sinusoid Benchmark Fatal Error Summary               \n');
            fprintf('==============================================================\n\n');
            failed_operations = readtable(testCase.failedOpFileName);
            disp(failed_operations);
            
        end
    end

end
