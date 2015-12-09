function [ord,R,keepers] = BF_ClusterReorder(dataMatrix,distanceMetric,linkageMethod)
% BF_ClusterReorder     A clustered reordering of the rows of an input data matrix.
%
% Computes a reordering of the rows of an input data matrix (under a given
% distance metric), placing similar rows close together in the output
% permutation, ord.
%
% Alternatively, you can input a distance matrix for distanceMetric if
% pre-computed.

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check Inputs:
%-------------------------------------------------------------------------------
if nargin < 2
    distanceMetric = 'corr'; % correlation distances by default
end
if nargin < 3
    linkageMethod = 'average'; % average linkage by default
end

% ------------------------------------------------------------------------------
%% Do linkage:
% ------------------------------------------------------------------------------
if ischar(distanceMetric)
    % Specify a distance metric as an input to pdist/BF_pdist
    R = BF_pdist(dataMatrix,distanceMetric);
else
    % Put the pre-computed distance matrix in the second input: distanceMetric
    R = distanceMetric;
end

% Convert pairwise distances to matrix using squareform if stored as a vector:
if size(R,1)==1 || size(R,2)==1
    R = squareform(R);
end

if any(isnan(R(:)))
    % Remove NaNs:
    [R,keepers] = BF_RemoveNaN_DistMat(R);
    fprintf(1,'***CAUTION: Removed %u bad features from the distance matrix\n', ...
                        sum(keepers==0));
else
    keepers = ones(length(R),1);
end

if size(R,1)==size(R,2)
    R = squareform(R); % Convert back to vector for linkage to work properly
end
links = linkage(R,linkageMethod);

% ------------------------------------------------------------------------------
%% Get the optimal dendrogram reordering:
% ------------------------------------------------------------------------------
figure('color','w');
set(gcf,'Visible','off'); % suppress figure output
if sqrt(length(R)) < 2000 % small enough to try optimalleaforder
    try
        ord = optimalleaforder(links,R); % NEW!
        [~,~,ord] = dendrogram(links,0,'r',ord);
        fprintf(1,'Using optimalleaforder reordering!\n');
    catch
        fprintf(1,'Using dendrogram reordering.\n');
        [~,~,ord] = dendrogram(links,0);
    end
else
    fprintf(1,'Too big for optimalleaforder, using dendrogram.\n');
    [~,~,ord] = dendrogram(links,0);
end
close; % close the invisible figure used for the dendrogram

if ~all(keepers==1)
    keepers = find(keepers);
    ord = keepers(ord); % convert to indicies of the input matrix
end

end
