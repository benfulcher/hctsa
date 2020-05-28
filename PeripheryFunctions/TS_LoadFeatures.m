function TS_LoadFeatures(whatDataFile,TS_DataMatToLoad,TS_CalcTimeToLoad,TS_QualityToLoad,INP_opsToLoad,INP_tsToLoad)

load(whatDataFile);

if isempty(INP_opsToLoad)
    myOpIds = 1:size(TS_DataMat,2);
else
    myOpIds = TS_GetIDs(INP_opsToLoad,whatDataFile,'ops','Name');
end
if isempty(INP_tsToLoad)
    myTSIds = 1:size(TS_DataMat,1);
else
    myTSIds = TS_GetIDs(INP_tsToLoad,whatDataFile,'ts','Name');
end

TS_DataMat(myTSIds,myOpIds) = TS_DataMatToLoad;
TS_CalcTime(myTSIds,myOpIds) = TS_CalcTimeToLoad;
TS_Quality(myTSIds,myOpIds) = TS_QualityToLoad;

save(whatDataFile,'TS_DataMat','TS_CalcTime','TS_Quality','-append');