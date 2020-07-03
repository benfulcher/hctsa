function TS_LoadFeatures(whatDataFile,TS_DataMatToLoad,TS_CalcTimeToLoad,TS_QualityToLoad,INP_opsToLoad,INP_tsToLoad)
% TS_LoadFeatures   Load some pre-computed features into an hctsa data file
%                     if using externally computed values, e.g., from Catch22.
%
%---INPUTS:
% whatDataFile, the input/output hctsa file (already created with TS_Init)
% TS_DataMatToLoad, the feature matrix to load, as a time series-by-features matrix
% TS_CalcTimeToLoad, the calculation time for each operation, same size as TS_DataMatToLoad
% TS_QualityToLoad, the quality for each operation (zero means no issues), same size as TS_DataMatToLoad
% INP_opsToLoad, the names of the hctsa operations, calculation and quality matrices to be loaded (see INP_ops.txt)
% INP_tsToLoad, the names of the time series that are in whatDataFile that we're loading
%
% N.B. the whatDataFile should contain all operations specified in TS_opsToLoad

load(whatDataFile);

% Which operations are we loading?
if isempty(INP_opsToLoad)
    myOpIDs = 1:size(TS_DataMat,2);
else
    myOpIDs = TS_GetIDs(INP_opsToLoad,whatDataFile,'ops','Name');
end

% Which time series are we loading?
if isempty(INP_tsToLoad)
    myTSIDs = 1:size(TS_DataMat,1);
else
    myTSIDs = TS_GetIDs(INP_tsToLoad,whatDataFile,'ts','Name');
end

TS_DataMat(myTSIDs,myOpIDs) = TS_DataMatToLoad;
TS_CalcTime(myTSIDs,myOpIDs) = TS_CalcTimeToLoad;
TS_Quality(myTSIDs,myOpIDs) = TS_QualityToLoad;

save(whatDataFile,'TS_DataMat','TS_CalcTime','TS_Quality','-append');

end
