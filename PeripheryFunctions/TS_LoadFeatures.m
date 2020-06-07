function TS_LoadFeatures(whatDataFile,TS_DataMatToLoad,TS_CalcTimeToLoad,TS_QualityToLoad,INP_opsToLoad,INP_tsToLoad)

load(whatDataFile);

if isempty(INP_opsToLoad)
    myOpIDs = 1:size(TS_DataMat,2);
else
    myOpIDs = TS_GetIDs(INP_opsToLoad,whatDataFile,'ops','Name');
end
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
