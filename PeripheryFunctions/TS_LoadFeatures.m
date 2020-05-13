function TS_LoadFeatures(TS_DataMat,TS_CalcTime,TS_Quality,whatDataFile)

save(whatDataFile,'TS_DataMat','TS_CalcTime','TS_Quality','-append');