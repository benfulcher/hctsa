function [T,XTrainNorm] = BF_FitTransform(X,spec,Operations,y,beVerbose)
% BF_FitTransform   Fits a feature filter and normalization on training data.
%
% Fits, on the rows of X only, the feature (column) steps of TS_Normalize:
% (i) features with too few good values are removed, (ii) constant features are
% removed, (iii) optionally, features that are constant within any class are
% removed, (iv) each remaining feature's normalization is fitted
% (BF_NormalizeMatrix), and (v) features that are all-NaN or constant after
% normalization are removed.
% The result is applied to any data (e.g., test rows or new data) by
% BF_ApplyTransform, which fits nothing.
%
% The row filter of TS_Normalize (removing time series with too few good
% values) uses only each time series' own values, so it is not part of the
% transformation: apply it separately.
%
%---USAGE:
% T = BF_FitTransform(XTrain,spec,Operations);
% XTestNorm = BF_ApplyTransform(XTest,T,Operations);
%
%---INPUTS:
% X, the data to fit on (training rows x features). Infs are treated as missing
%       values (NaN), as in TS_Normalize.
% spec, the normalization specification, a structure with (optional) fields:
%   normFunction: a normalization method of BF_NormalizeMatrix, or 'none' for no
%                 normalization. Default: 'mixedSigmoid'.
%   filterOptions: [row proportion, column proportion] of good values required
%                 (as in TS_Normalize). Only the column proportion is used here.
%                 Default: [0.70,1].
%   classVarFilter: whether to remove features with near-zero variance within any
%                 class (requires y). Default: false.
%   testNaN: how BF_ApplyTransform treats missing values of kept features:
%                 'keep' (they stay NaN), or 'impute' (they are replaced by the
%                 median of the normalized training values of that feature).
%                 Default: 'keep'.
%   (spec can be empty, or a normalizationInfo structure saved by TS_Normalize)
% Operations, a table with the ID and Name of each column of X, or a vector of
%       column IDs. Default: IDs 1, 2, ..., and empty names.
% y, the class labels of the rows of X (only used if spec.classVarFilter).
% beVerbose, whether to describe each step. Default: false.
%
%---OUTPUTS:
% T, the fitted transformation, a structure with fields:
%   spec: the specification (with defaults filled in)
%   numTrain: the number of rows it was fitted on
%   featureParams: a table with a row for each feature kept (in the column order
%       of X), with variables ID and Name, and (unless spec.normFunction is
%       'none') the parameters of its normalization: transform, location, scale,
%       rescaleMin, rescaleMax (cf. the normParams output of BF_NormalizeMatrix).
%   trainMedian: the median of the normalized training values of each feature
%       kept (a row vector, used when spec.testNaN is 'impute').
% XTrainNorm, the features of X that are kept, transformed (the same as
%       BF_ApplyTransform(X,T,Operations,false)).

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
if nargin < 2
    spec = [];
end
spec = FillSpec(spec);
[numTrain,numFeatures] = size(X);
if nargin < 3 || isempty(Operations)
    Operations = table((1:numFeatures)',repmat({''},numFeatures,1),'VariableNames',{'ID','Name'});
elseif isnumeric(Operations)
    Operations = table(Operations(:),repmat({''},numFeatures,1),'VariableNames',{'ID','Name'});
end
if height(Operations)~=numFeatures
    error('Operations has %u rows, but X has %u columns',height(Operations),numFeatures);
end
if nargin < 4
    y = [];
end
if nargin < 5
    beVerbose = false;
end
if numTrain < 2
    error('At least two rows are needed to fit the transformation (%u given)',numTrain);
end

X(~isfinite(X)) = NaN; % all special values are treated as missing
keepCol = (1:numFeatures)'; % indices of the columns of X that are kept

%-------------------------------------------------------------------------------
%% (i) Features with too few good values
%-------------------------------------------------------------------------------
colThreshold = spec.filterOptions(2);
if colThreshold > 0
    isGood = (1 - mean(isnan(X),1) >= colThreshold)';
    if ~any(isGood)
        error('No operations had more than %4.2f%% good values.\nSet a more lenient threshold.',...
                            colThreshold*100)
    end
    if beVerbose
        if all(isGood)
            fprintf(1,'All %u operations have greater than %4.2f%% good values. Keeping them all.\n',...
                            numFeatures,colThreshold*100);
        else
            fprintf(1,'Removing %u operations with fewer than %4.2f%% good values: from %u to %u.\n',...
                            sum(~isGood),colThreshold*100,numFeatures,sum(isGood));
        end
    end
    keepCol = keepCol(isGood);
end

%-------------------------------------------------------------------------------
%% (ii) Constant features
%-------------------------------------------------------------------------------
isConstant = (std(X(:,keepCol),0,'omitnan') < 10*eps)';
if all(isConstant)
    error('All %u operations produced constant outputs on the %u time series?!',...
                        length(isConstant),numTrain)
elseif any(isConstant)
    if beVerbose
        fprintf(1,'Removed %u operations with near-constant outputs: from %u to %u.\n',...
                        sum(isConstant),length(isConstant),sum(~isConstant));
    end
    keepCol = keepCol(~isConstant);
elseif beVerbose
    fprintf(1,'No operations had near-constant outputs on the dataset\n');
end

%-------------------------------------------------------------------------------
%% (iii) Features with near-zero variance within a class
%-------------------------------------------------------------------------------
if spec.classVarFilter
    if isempty(y)
        error('Class labels are required to filter on class variance');
    end
    y = y(:);
    if ~iscategorical(y)
        y = categorical(y);
    end
    classNames = categories(y);
    numClasses = length(classNames);
    classVars = zeros(numClasses,length(keepCol));
    for i = 1:numClasses
        classVars(i,:) = std(X(y==classNames{i},keepCol),0,'omitnan');
    end
    zeroClassVar = any(classVars < 10*eps,1)';
    if all(zeroClassVar)
        error('All %u operations produced near-constant class-wise outputs?!',length(zeroClassVar))
    elseif any(zeroClassVar)
        if beVerbose
            fprintf(1,'Removed %u operations with near-constant class-wise outputs: from %u to %u.\n',...
                        sum(zeroClassVar),length(zeroClassVar),sum(~zeroClassVar));
        end
        keepCol = keepCol(~zeroClassVar);
    end
end

%-------------------------------------------------------------------------------
%% (iv) Fit the normalization of each feature
%-------------------------------------------------------------------------------
XTrainNorm = X(:,keepCol);
if beVerbose
    numBad = sum(isnan(XTrainNorm(:)));
    if numBad==0
        fprintf(1,'\n(post-filtering): No special-valued entries in the %ux%u data matrix!\n',...
                    size(XTrainNorm,1),size(XTrainNorm,2));
    else
        fprintf(1,['\n(post-filtering): %u special-valued entries (%4.2f%%) remain in the' ...
                    ' %ux%u data matrix.\n'],numBad,numBad/numel(XTrainNorm)*100,...
                    size(XTrainNorm,1),size(XTrainNorm,2));
    end
end
doNormalize = ~ismember(spec.normFunction,{'none','nothing'});
if doNormalize
    if beVerbose
        fprintf(1,'Normalizing a %u x %u object. Please be patient...\n',size(XTrainNorm,1),size(XTrainNorm,2));
    end
    [XTrainNorm,normParams] = BF_NormalizeMatrix(XTrainNorm,spec.normFunction);
    featureParams = normParams.featureParams;
elseif beVerbose
    fprintf(1,'You specified ''%s'', so NO NORMALIZING IS ACTUALLY BEING DONE!!!\n',spec.normFunction);
end

%-------------------------------------------------------------------------------
%% (v) Features that are all-NaN or constant after normalization
%-------------------------------------------------------------------------------
isAllNaN = (mean(isnan(XTrainNorm),1)==1)';
if all(isAllNaN)
    error('After normalization, all columns were bad-values... :(');
elseif any(isAllNaN) && beVerbose
    fprintf(1,'We just removed %u all-NaN columns introduced from %s normalization.\n',...
                    sum(isAllNaN),spec.normFunction);
end
XTrainNorm = XTrainNorm(:,~isAllNaN);
keepCol = keepCol(~isAllNaN);
if doNormalize
    featureParams = featureParams(~isAllNaN,:);
end

isConstant = (std(XTrainNorm,0,'omitnan') < 10*eps)';
if any(isConstant)
    if beVerbose
        fprintf(1,'%u operations had near-constant outputs after filtering: from %u to %u.\n',...
                    sum(isConstant),length(isConstant),sum(~isConstant));
    end
    XTrainNorm = XTrainNorm(:,~isConstant);
    keepCol = keepCol(~isConstant);
    if doNormalize
        featureParams = featureParams(~isConstant,:);
    end
end

%-------------------------------------------------------------------------------
%% Store the transformation
%-------------------------------------------------------------------------------
T = struct();
T.spec = spec;
T.numTrain = numTrain;
if doNormalize
    T.featureParams = [Operations(keepCol,{'ID','Name'}),featureParams];
else
    T.featureParams = Operations(keepCol,{'ID','Name'});
end
T.trainMedian = median(XTrainNorm,1,'omitnan');

end

%-------------------------------------------------------------------------------
function spec = FillSpec(spec)
    % Normalization specification with defaults for unset fields
    if isempty(spec)
        spec = struct();
    end
    if ~isfield(spec,'normFunction') || isempty(spec.normFunction)
        spec.normFunction = 'mixedSigmoid';
    end
    if ~isfield(spec,'filterOptions') || isempty(spec.filterOptions)
        spec.filterOptions = [0.70,1];
    end
    if ~isfield(spec,'classVarFilter') || isempty(spec.classVarFilter)
        spec.classVarFilter = false;
    end
    if ~isfield(spec,'testNaN') || isempty(spec.testNaN)
        spec.testNaN = 'keep';
    end
    if ~ismember(spec.testNaN,{'keep','impute'})
        error('spec.testNaN must be ''keep'' or ''impute''');
    end
    % Keep only the fields of the specification (e.g., from a normalizationInfo structure):
    spec = struct('normFunction',spec.normFunction,'filterOptions',spec.filterOptions,...
                    'classVarFilter',logical(spec.classVarFilter),'testNaN',spec.testNaN);
end
