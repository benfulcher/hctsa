function TS_LoadFeatures(whatDataFile,TS_DataMatToLoad,TS_CalcTimeToLoad,TS_QualityToLoad,INP_opsToLoad,INP_tsToLoad)
% TS_LoadFeatures   Load some pre-computed features into an hctsa data file
%                     if using externally computed values, e.g., from Catch22.
%
%---INPUTS:
% whatDataFile, the input/output hctsa file (already created with TS_Init)
% TS_DataMatToLoad, the feature matrix to load, as a time series-by-features
%               matrix
% TS_CalcTimeToLoad, the calculation time for each operation, same size as
%                       TS_DataMatToLoad
% TS_QualityToLoad, the quality for each operation (zero means no issues), same
%                       size as TS_DataMatToLoad
% INP_opsToLoad, the names of the hctsa operations, calculation and quality
%                   matrices to be loaded (see INP_ops.txt)
% INP_tsToLoad, the names of the time series that are in whatDataFile that
%                   we're loading
%
% whatDataFile should contain all operations specified in TS_opsToLoad

% ------------------------------------------------------------------------------
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

% whatDataFile must be a filename

% Load in the data:
TS_DataMat = TS_GetFromData(whatDataFile,'TS_DataMat');
TS_CalcTime = TS_GetFromData(whatDataFile,'TS_CalcTime');
TS_Quality = TS_GetFromData(whatDataFile,'TS_Quality');

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

% Fill the data:
TS_DataMat(myTSIDs,myOpIDs) = TS_DataMatToLoad;
TS_CalcTime(myTSIDs,myOpIDs) = TS_CalcTimeToLoad;
TS_Quality(myTSIDs,myOpIDs) = TS_QualityToLoad;

% Save to the data file
save(whatDataFile,'TS_DataMat','TS_CalcTime','TS_Quality','-append');

end
