function TS_CompareSpeed(dataInput)
% Compares the speed of calculation between different features

if isnumeric(dataInput)
    [featureVector,calcTimes,calcQuality,Operations,MasterOperations] = TS_CalculateFeatureVector(dataInput);
end

% sort
[~,ix] = sort(calcTimes,'descend');
notNaN = ~isnan(calcTimes(ix));
ix = ix(notNaN);

keyboard

for i = 1:100
    ixi = ix(i);
    fprintf(1,'[%u] %s, %s: %s\n',Operations(ixi).ID,Operations(ixi).Name,Operations(ixi).CodeString,BF_thetime(calcTimes(ixi)));
end

end
