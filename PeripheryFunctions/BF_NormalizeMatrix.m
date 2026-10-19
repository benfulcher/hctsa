function [dataMatrixNorm,normParams] = BF_NormalizeMatrix(dataMatrix,normMethod,isTraining)
% BF_NormalizeMatrix    Normalizes all columns of an input matrix.
%
% Each column is transformed separately. The parameters of each transformation
% are fitted on the training rows and the transformation is then applied to all
% rows. The fitted parameters can be returned and applied to new data later.
%
%---USAGE:
% XNorm = BF_NormalizeMatrix(X,'mixedSigmoid'); % fit on all rows, transform all rows
% XNorm = BF_NormalizeMatrix(X,'zscore',isTraining); % fit on training rows only
% [XNorm,normParams] = BF_NormalizeMatrix(X,'mixedSigmoid'); % also return the parameters
% XNewNorm = BF_NormalizeMatrix(XNew,normParams); % apply fitted parameters to new data
%
%---INPUTS:
% dataMatrix, the input data matrix (observations as rows, features as columns)
% normMethod, the normalization method to use:
%   'subtractMean': x - mean(x)
%   'maxmin': linear rescaling to the unit interval, (x - min(x))/(max(x) - min(x))
%   'zscore': (x - mean(x))/std(x)
%   'sigmoid': 1/(1 + exp(-(x - mean(x))/std(x)))
%   'scaledSigmoid': 'sigmoid', then a linear rescaling to the unit interval
%   'robustSigmoid': outlier-robust sigmoid, 1/(1 + exp(-(x - median(x))/(iqr(x)/1.35)))
%   'scaledRobustSigmoid': 'robustSigmoid', then a linear rescaling to the unit interval
%   'mixedSigmoid': 'scaledRobustSigmoid', or 'scaledSigmoid' for a column with iqr = 0
%   'scaledsigmoid5q': 'scaledSigmoid' with the mean and std computed from values
%                       between the 5th and 95th percentiles only
%   OR the normParams output of an earlier call: its fitted parameters are then
%   applied to the columns of dataMatrix (nothing is fitted; isTraining is not used).
%   Default: 'sigmoid'.
% isTraining, learn the normalization parameters just on this training portion
%             of data, then apply it to the full dataset (required for training
%             /testing procedures where the testing data has to remain unseen).
%             Should be input as a logical of the same size as the number of rows
%             in the dataMatrix. Default: all rows.
%
%---OUTPUTS:
% dataMatrixNorm, the normalized matrix
% normParams, the fitted parameters, a structure with fields:
%   normMethod, the normalization method
%   featureParams, a table with a row for each column of dataMatrix, with variables:
%       transform: the transformation applied to that column: 'subtractMean',
%                  'maxmin', 'zscore', 'sigmoid', 'scaledSigmoid', 'robustSigmoid',
%                  'scaledRobustSigmoid', 'constant' (values set to 0: constant
%                  training values), or 'allNaN' (values set to NaN: no training
%                  values, or a scale of zero for a robust sigmoid or 'scaledsigmoid5q').
%       location: the mean (median for robust sigmoids, minimum for 'maxmin')
%       scale: the std (iqr for robust sigmoids, max - min for 'maxmin')
%       rescaleMin, rescaleMax: for the scaled sigmoids, the range of the sigmoid
%                  output on the training rows, used to rescale to the unit interval
%                  (all values are set to 0 if they are equal; no rescaling is done
%                  if they are NaN, i.e., no training output is a real number).
%
% NaNs are ignored when fitting (only real data is used for the normalization),
% and stay NaN in the output. Infs are not treated specially (TS_Normalize
% converts them to NaN before normalizing).

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
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%% Check inputs, set defaults
% ------------------------------------------------------------------------------

if nargin < 2 || isempty(normMethod)
    fprintf(1,'We''re normalizing using sigmoid transform by default\n');
    normMethod = 'sigmoid';
end

numObs = size(dataMatrix,1);
numFeatures = size(dataMatrix,2);

% ------------------------------------------------------------------------------
%% Apply previously fitted parameters
% ------------------------------------------------------------------------------
if isstruct(normMethod)
    normParams = normMethod;
    if ~isfield(normParams,'featureParams') || height(normParams.featureParams)~=numFeatures
        error(['Normalization parameters must be the normParams output of BF_NormalizeMatrix,' ...
                    ' with one row of featureParams for each of the %u columns'],numFeatures);
    end
    fp = normParams.featureParams;
    dataMatrixNorm = zeros(size(dataMatrix));
    for i = 1:numFeatures
        dataMatrixNorm(:,i) = ApplyTransform(dataMatrix(:,i),fp.transform{i},fp.location(i),...
                                fp.scale(i),fp.rescaleMin(i),fp.rescaleMax(i));
    end
    return
end

if nargin < 3 || isempty(isTraining)
    isTraining = true(numObs,1);
end
isTraining = isTraining(:); % should be a column vector
if ~islogical(isTraining) || length(isTraining)~=numObs
    error('isTraining should be a logical vector with one entry for each of the %u rows',numObs);
end

% ------------------------------------------------------------------------------
%% Fit the parameters of each column on the training rows, and transform
% ------------------------------------------------------------------------------
transform = cell(numFeatures,1);
location = nan(numFeatures,1);
scale = nan(numFeatures,1);
rescaleMin = nan(numFeatures,1);
rescaleMax = nan(numFeatures,1);
dataMatrixNorm = zeros(size(dataMatrix));

for i = 1:numFeatures
    x = dataMatrix(:,i);
    xTrain = x(~isnan(x) & isTraining); % only learn from real training data

    switch normMethod
        case 'subtractMean'
            transform{i} = 'subtractMean';
            location(i) = mean(xTrain);

        case 'maxmin'
            % Linear rescaling to the unit interval
            if isempty(xTrain)
                transform{i} = 'allNaN';
            elseif min(xTrain)==max(xTrain)
                % There is no variation -- set all to zero -- attempts to rescale will blow up
                transform{i} = 'constant';
            else
                transform{i} = 'maxmin';
                location(i) = min(xTrain);
                scale(i) = max(xTrain) - min(xTrain);
            end

        case 'zscore'
            transform{i} = 'zscore';
            location(i) = mean(xTrain);
            scale(i) = std(xTrain);

        case {'sigmoid','scaledSigmoid'}
            % Standard sigmoidal transformation (optionally rescaled to the unit interval)
            transform{i} = normMethod;
            location(i) = mean(xTrain);
            scale(i) = std(xTrain);

        case {'robustSigmoid','scaledRobustSigmoid'}
            % An outlier-robust sigmoid (optionally rescaled to the unit interval)
            % It is not defined if iqr = 0
            location(i) = median(xTrain);
            scale(i) = iqr(xTrain);
            if scale(i)==0
                transform{i} = 'allNaN';
            else
                transform{i} = normMethod;
            end

        case 'mixedSigmoid'
            % Uses a scaled sigmoid if iqr=0; a scaled, outlier-robust sigmoid otherwise
            if isempty(xTrain)
                % Everything a NaN, kept at NaN:
                transform{i} = 'allNaN';
            elseif max(xTrain)==min(xTrain)
                % A constant column is set to 0:
                transform{i} = 'constant';
            elseif iqr(xTrain)==0
                % iqr of data is zero: perform a normal sigmoidal transformation:
                transform{i} = 'scaledSigmoid';
                location(i) = mean(xTrain);
                scale(i) = std(xTrain);
            else
                % Perform an outlier-robust version of the sigmoid:
                transform{i} = 'scaledRobustSigmoid';
                location(i) = median(xTrain);
                scale(i) = iqr(xTrain);
            end

        case 'scaledsigmoid5q'
            % Scaled sigmoid, with mean and std computed from values within the
            % 5th--95th quantile range only
            qs = quantile(xTrain,[0.05,0.95]);
            qr = (xTrain>=qs(1) & xTrain<=qs(2)); % quantile range
            location(i) = mean(xTrain(qr));
            scale(i) = std(xTrain(qr));
            if isempty(xTrain) || scale(i)==0
                transform{i} = 'allNaN'; % avoid +/- Infs
            else
                transform{i} = 'scaledSigmoid';
            end

        otherwise
            error('Invalid normalization method ''%s''', normMethod)
    end

    % Range of the sigmoid output on the training rows, for rescaling to the unit interval:
    if ismember(transform{i},{'scaledSigmoid','scaledRobustSigmoid'})
        xhat = ApplyTransform(x,transform{i},location(i),scale(i),NaN,NaN);
        xhatTrain = xhat(~isnan(xhat) & isTraining);
        if ~isempty(xhatTrain)
            rescaleMin(i) = min(xhatTrain);
            rescaleMax(i) = max(xhatTrain);
        end
    end

    dataMatrixNorm(:,i) = ApplyTransform(x,transform{i},location(i),scale(i),...
                                            rescaleMin(i),rescaleMax(i));
end

if nargout > 1
    normParams = struct('normMethod',normMethod,'featureParams',...
                    table(transform,location,scale,rescaleMin,rescaleMax));
end

%-------------------------------------------------------------------------------
% Normalization functions:
%-------------------------------------------------------------------------------
function xhat = ApplyTransform(x,transform,location,scale,rescaleMin,rescaleMax)
    % Apply the transformation of a single column, given its fitted parameters
    switch transform
        case 'subtractMean'
            xhat = x - location;
        case {'maxmin','zscore'}
            xhat = (x - location)/scale;
        case {'sigmoid','scaledSigmoid'}
            % Classic sigmoidal transformation:
            xhat = 1./(1 + exp(-(x - location)/scale));
        case {'robustSigmoid','scaledRobustSigmoid'}
            % Outlier-robust sigmoid:
            xhat = 1./(1 + exp(-(x - location)/(scale/1.35)));
        case 'constant'
            xhat = zeros(size(x));
            xhat(isnan(x)) = NaN;
        case 'allNaN'
            xhat = nan(size(x));
        otherwise
            error('Unknown transformation ''%s''',transform);
    end

    % Rescale the sigmoid output to the unit interval (using its training range):
    if ismember(transform,{'scaledSigmoid','scaledRobustSigmoid'}) && ~isnan(rescaleMin)
        if rescaleMin==rescaleMax
            % There is no variation -- set all to zero -- attempts to rescale will blow up
            xhat(~isnan(xhat)) = 0;
        else
            xhat = (xhat - rescaleMin)/(rescaleMax - rescaleMin);
        end
    end
end

end
