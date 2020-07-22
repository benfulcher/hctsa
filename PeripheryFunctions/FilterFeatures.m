function [TS_DataMat,Operations] = FilterFeatures(TS_DataMat,Operations,cfnParams)
% FilterFeatures filters down to a reduced feature set
%-------------------------------------------------------------------------------
if isempty(cfnParams.reducedFeatureSet) || strcmp(cfnParams.reducedFeatureSet,'all')
    return
else
    opIDs = GiveMeFeatureSet(cfnParams.reducedFeatureSet,Operations);
    keepMe = ismember(Operations.ID,opIDs);
    Operations = Operations(keepMe,:);
    TS_DataMat = TS_DataMat(:,keepMe);
    fprintf(1,'Only using %u features from the %s set...\n',sum(keepMe),cfnParams.reducedFeatureSet);
end

% Check there are some left:
numFeatures = height(Operations);
if numFeatures==0
    error('No features to do classification with...');
end

end
