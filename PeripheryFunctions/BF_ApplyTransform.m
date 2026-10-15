function [XNorm,isImputed] = BF_ApplyTransform(X,T,Operations,doImpute)
% BF_ApplyTransform   Applies a transformation fitted by BF_FitTransform.
%
% Selects the features kept by the transformation (by ID) and normalizes each
% with its stored parameters. Nothing is fitted to X, so test data (or new data)
% are transformed exactly as the training data were.
%
% Values outside the range of the training data are transformed with the same
% parameters: sigmoid outputs stay in the unit interval, but 'maxmin' and the
% rescaling of the scaled sigmoids (based on the training range) can give values
% outside [0,1] (they are not clipped).
%
%---USAGE:
% T = BF_FitTransform(XTrain,spec,Operations);
% XTestNorm = BF_ApplyTransform(XTest,T,Operations);
%
%---INPUTS:
% X, the data (rows x features). Infs are treated as missing values (NaN).
% T, the transformation (output of BF_FitTransform).
% Operations, a table with the ID of each column of X (or a vector of IDs).
%       Default (empty): the columns of X are the features of T, in order.
% doImpute, whether to replace missing values of the kept features by the median
%       of the normalized training values of that feature.
%       Default: true if T.spec.testNaN is 'impute', false otherwise.
%
%---OUTPUTS:
% XNorm, the transformed data: one column for each feature of T (in the order of
%       T.featureParams).
% isImputed, a logical matrix (the size of XNorm) marking the imputed values.

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
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

%-------------------------------------------------------------------------------
%% Check inputs and set defaults
%-------------------------------------------------------------------------------
featureParams = T.featureParams;
numFeatures = height(featureParams);
if nargin < 4 || isempty(doImpute)
    doImpute = strcmp(T.spec.testNaN,'impute');
end

%-------------------------------------------------------------------------------
%% Select the features of the transformation
%-------------------------------------------------------------------------------
if nargin < 3 || isempty(Operations)
    if size(X,2)~=numFeatures
        error('X has %u columns, but the transformation has %u features',size(X,2),numFeatures);
    end
else
    if istable(Operations)
        colIDs = Operations.ID;
    else
        colIDs = Operations(:);
    end
    if length(colIDs)~=size(X,2)
        error('%u column IDs given for a matrix with %u columns',length(colIDs),size(X,2));
    end
    [isFound,colInd] = ismember(featureParams.ID,colIDs);
    if ~all(isFound)
        error('%u features of the transformation are not in the data (e.g., ID %u)',...
                    sum(~isFound),featureParams.ID(find(~isFound,1)));
    end
    X = X(:,colInd);
end
X(~isfinite(X)) = NaN; % all special values are treated as missing

%-------------------------------------------------------------------------------
%% Normalize with the stored parameters
%-------------------------------------------------------------------------------
if ismember('transform',featureParams.Properties.VariableNames)
    normParams = struct('normMethod',T.spec.normFunction,'featureParams',featureParams);
    XNorm = BF_NormalizeMatrix(X,normParams);
else
    % No normalization ('none'):
    XNorm = X;
end

%-------------------------------------------------------------------------------
%% Impute missing values
%-------------------------------------------------------------------------------
if doImpute
    medianMat = repmat(T.trainMedian,size(XNorm,1),1);
    isImputed = isnan(XNorm) & ~isnan(medianMat);
    XNorm(isImputed) = medianMat(isImputed);
else
    isImputed = false(size(XNorm));
end

end
