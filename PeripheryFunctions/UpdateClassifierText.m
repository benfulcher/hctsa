function cfnParams = UpdateClassifierText(cfnParams)

switch cfnParams.whatClassifier
    case {'linear','linclass','fast_linear'}
        cfnParams.classifierText = 'linear classifier';
    case 'diaglinear'
        cfnParams.classifierText = 'linear naive bayes classifier';
    case {'svm','svm_linear'}
        cfnParams.classifierText = 'linear SVM classifier';
    otherwise
        cfnParams.classifierText = cfnParams.whatClassifier;
end

end
