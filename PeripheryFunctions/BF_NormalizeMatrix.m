function dataMatrixNorm = BF_NormalizeMatrix(dataMatrix,normMethod,isTraining)
% BF_NormalizeMatrix    Normalizes all columns of an input matrix.
%
%---INPUTS:
% dataMatrix, the input data matrix
% normMethod, the normalization method to use (see body of the code for options)
% isTraining, learn the normalization parameters just on this training portion
%             of data, then apply it to the full dataset (required for training
%             /testing procedures where the testing data has to remain unseen).
%             Should be input as a logical of the same size as the number of rows
%             in the dataMatrix.
%
%---OUTPUT:
% dataMatrixNorm, the normalized matrix
%
% NaNs are ignored -- only real data is used for the normalization

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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

if nargin < 3
    isTraining = true(size(dataMatrix,1),1);
end
if size(isTraining,2) > size(isTraining,1)
    isTraining = isTraining'; % should be a column vector
end
if size(isTraining,1)~=1
    error('Logical specifying training indices should be a vector');
end
%-------------------------------------------------------------------------------

numFeatures = size(dataMatrix,2);

% ------------------------------------------------------------------------------
% Normalize each column according to the specified normalizing transformation
% ------------------------------------------------------------------------------

dataMatrixNorm = zeros(size(dataMatrix));

switch normMethod
    case 'subtractMean'
        % Subtract the mean:
        dataMatrixNorm = bsxfun(@minus,dataMatrix,mean(dataMatrix));

    case 'maxmin'
        % Linear rescaling to the unit interval
        for i = 1:numFeatures % cycle through the operations
            dataMatrixNorm(:,i) = UnityRescale(dataMatrix(:,i));
        end

    case 'zscore'
        for i = 1:numFeatures
            dataMatrixNorm(:,i) = ZScore(dataMatrix(:,i));
        end

    case 'robustSigmoid'
        % A outlier-robust sigmoid
        for i = 1:numFeatures % cycle through the features
            dataMatrixNorm(:,i) = RobustSigmoid(dataMatrix(:,i),0);
        end

    case 'scaledRobustSigmoid'
        % A scaled, outlier-robust sigmoid
        % Problem is that if the iqr=0, this is not defined

        for i = 1:numFeatures % cycle through the features
            dataMatrixNorm(:,i) = RobustSigmoid(dataMatrix(:,i),1);
        end

    case 'sigmoid'
        % Standard sigmoidal transformation
        for i = 1:numFeatures % cycle through the features
            dataMatrixNorm(:,i) = Sigmoid(dataMatrix(:,i),0);
        end

    case 'scaledSigmoid'
        % Standard sigmoid transform, then a rescaling to the unit interval
        dataMatrixNorm = zeros(size(dataMatrix));
        for i = 1:numFeatures % cycle through the features
            dataMatrixNorm(:,i) = Sigmoid(dataMatrix(:,i),1);
        end

    case 'mixedSigmoid'
        % Uses a scaled sigmoid if iqr=0; a scaled, outlier-robust sigmoid otherwise
        % Uses only non-NaNs

        % Outlier-adjusted sigmoid:
        for i = 1:numFeatures % cycle through columns
            theColumn = dataMatrix(:,i);
            if max(dataMatrix(:,i))==min(theColumn)
                % A constant column is set to 0:
                dataMatrixNorm(:,i) = 0;
            elseif all(isnan(theColumn))
                % Everything a NaN, kept at NaN:
                dataMatrixNorm(:,i) = NaN;
            elseif iqr(dataMatrix(~isnan(theColumn),i))==0
                % iqr of data is zero: perform a normal sigmoidal transformation:
                dataMatrixNorm(:,i) = Sigmoid(dataMatrix(:,i),1);
            else
                % Perform an outlier-robust version of the sigmoid:
                dataMatrixNorm(:,i) = RobustSigmoid(dataMatrix(:,i),1);
            end
        end

    case 'scaledsigmoid5q'
        % First caps at 5th and 95th quantile, then does scaled sigmoid
        for i = 1:numFeatures % cycle through the features
            goodRows = ~isnan(dataMatrix(:,i));
            FF = dataMatrix(goodRows,i);
            qs = quantile(FF,[0.05,0.95]);
            qr = (FF>=qs(1) & FF<=qs(2)); % quantile range
            % calculate mean and std based on quantile range only
            meanF = mean(FF(qr));
            stdF = std(FF(qr));
            if stdF==0
                dataMatrixNorm(goodRows,i) = NaN; % avoid +/- Infs
            else
%                 kk = 1./(1+exp(-zscore(FF)));
                kk = 1./(1+exp(-(FF-meanF)/stdF));
                dataMatrixNorm(goodRows,i) = (kk-min(kk))/(max(kk)-min(kk));
            end
        end

    otherwise
        error('Invalid normalization method ''%s''', normMethod)
end


%-------------------------------------------------------------------------------
% Normalization functions:
%-------------------------------------------------------------------------------

function xhat = Sigmoid(x,doScale)
    % Classic sigmoidal transformation (optionally scaled to the unit interval)
    if nargin < 2, doScale = 1; end

    goodVals = (~isnan(x) & isTraining);
    meanX = mean(x(goodVals));
    stdX = std(x(goodVals));

    % Sigmoidal transformation:
    xhat = 1./(1 + exp(-(x-meanX)/stdX));

    % Rescale to unit interval:
    if doScale
        xhat = UnityRescale(xhat);
    end
end

function xhat = RobustSigmoid(x,doScale)
    % Outlier-adjusted sigmoid (optionally scaled to unit interval)
    if nargin < 2, doScale = 1; end

    goodVals = (~isnan(x) & isTraining);
    medianX = median(x(goodVals));
    iqrX = iqr(x(goodVals));

    if iqrX==0 % Can't apply an outlier-robust sigmoid meaningfully
        xhat = ones(size(x))*NaN;
    else
        % Outlier-robust sigmoid:
        xhat = 1./(1 + exp(-(x-medianX)/(iqrX/1.35)));
        if doScale % Rescale to the unit interval:
            xhat = UnityRescale(xhat);
        end
    end
end

function xhat = ZScore(x)
    goodVals = (~isnan(x) & isTraining);
    meanX = mean(x(goodVals));
    stdX = std(x(goodVals));
    xhat = (x-meanX)/stdX;
end

function xhat = UnityRescale(x)
    % Linearly rescale a data vector to unit interval:
    goodVals = (~isnan(x) & isTraining); % only work with non-NaN data
    if ~any(goodVals) % all NaNs:
        xhat = x; return
    end
    minX = min(x(goodVals));
    maxX = max(x(goodVals));
    if minX==maxX
        % There is no variation -- set all to zero -- attempts to rescale will blow up
        xhat = x;
        xhat(goodVals) = 0;
    else
        % There is some variation -- rescale to unit interval
        xhat = (x-minX)/(maxX-minX);
    end
end

end
