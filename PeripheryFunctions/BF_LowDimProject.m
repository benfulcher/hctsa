function [lowDimComponents,componentLabels] = BF_LowDimProject(dataMatrix,whatAlgorithm,numComponents,Operations)
% BF_LowDimProject Compute low-dimensional components of a data matrix.
%
%---INPUTS:
% dataMatrix, matrix of observations x variables
% whatAlgorithm, the dimensionality reduction algorithm ('pca','tnse')
% numComponents, the number of components to use
% Operations, table detailing the columns of dataMatrix (optional)

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

if nargin < 2
    whatAlgorithm = 'pca';
end
if nargin < 3
    numComponents = 2;
end
if nargin < 4
    Operations = [];
end
%-------------------------------------------------------------------------------
numFeatures = size(dataMatrix,2);
switch whatAlgorithm
case {'pca','PCA'}
    fprintf(1,'Calculating 2-dimensional principal components of the %u x %u data matrix...\n', ...
                        size(dataMatrix,1),size(dataMatrix,2));

    % Use pca to compute the first two principal components:
    % (project data into space of PC scores, Y)
    if ~any(isnan(dataMatrix))
        [pcCoeff,lowDimComponents,~,~,percVar] = pca(zscore(dataMatrix),'NumComponents',numComponents);
    else
        warning(sprintf(['Data matrix contains %.2g%% NaNs. Estimating covariances on remaining data...\n' ...
                    '(Could take some time...)'],100*mean(isnan(dataMatrix(:)))))
        % Data matrix contains NaNs; try the pairwise rows approximation to the
        % covariance matrix:
        [pcCoeff,lowDimComponents,~,~,percVar] = pca(BF_NormalizeMatrix(dataMatrix,'zscore'),...
                                    'Rows','pairwise','NumComponents',numComponents);
        % If this fails (covariance matrix not positive definite), can try the
        % (...,'algorithm','als') option in pca... (or toolbox for probabilistic PCA)
    end
    fprintf(1,'---Done.\n');

    %-------------------------------------------------------------------------------
    % Display the features loading strongly into the first two components:
    if ~isempty(Operations)
        numTopLoadFeat = min(numFeatures,20); % display this many features loading onto each PC
        LowDimDisplayTopLoadings(numTopLoadFeat,2,pcCoeff,lowDimComponents,dataMatrix,Operations);
    end

    % Text labels for each component:
    componentLabels = cell(2,1);
    for i = 1:2
        componentLabels{i} = sprintf('PC-%u (%.2f%% var)',i,percVar(i));
    end

case {'tSNE','tsne'}
    defaultNumPCs = 100;
    numPCAComponents = min(size(dataMatrix,2),defaultNumPCs);
    rng('default') % for reproducibility

    if numPCAComponents < size(dataMatrix,2)
        fprintf(1,['Computing a two-dimensional t-SNE embedding (using barnes-hut',...
                        ' approximation after %u-dim PC reduction) of the %u x %u data matrix...\n'], ...
                            numPCAComponents,size(dataMatrix,1),size(dataMatrix,2));
        lowDimComponents = tsne(BF_NormalizeMatrix(dataMatrix,'zscore'),'Algorithm','barneshut',...
                'Distance','euclidean','NumPCAComponents',numPCAComponents,'NumDimensions',numComponents);

    else
        fprintf(1,['Computing a two-dimensional t-SNE embedding (using barnes-hut',...
                        ' approximation) of the %u x %u data matrix...\n'], ...
                            size(dataMatrix,1),size(dataMatrix,2));
        lowDimComponents = tsne(BF_NormalizeMatrix(dataMatrix,'zscore'),'Algorithm','barneshut',...
                        'Distance','euclidean','NumDimensions',numComponents);
    end
    fprintf(1,'---Done.\n');
    % Axis labels for the plot:
    componentLabels = arrayfun(@(x)sprintf('tSNE-%u',x),1:numComponents,'UniformOutput',false);
otherwise
    error('Unknown dimensionality-reduction algorithm: %s',whatAlgorithm);
end

end
