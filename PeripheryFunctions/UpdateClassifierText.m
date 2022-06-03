function cfnParams = UpdateClassifierText(cfnParams)

switch cfnParams.whatClassifier
    case {'linear','linclass','fast-linear'}
        cfnParams.classifierText = 'simple linear classifier';
    case 'diaglinear'
        cfnParams.classifierText = 'linear naive Bayes classifier';
    case {'svm','svm-linear'}
        cfnParams.classifierText = 'linear SVM classifier';
    case {'svm-rbf'}
        cfnParams.classifierText = 'radial basis function SVM classifier';
    otherwise
        cfnParams.classifierText = cfnParams.whatClassifier;
end

end
