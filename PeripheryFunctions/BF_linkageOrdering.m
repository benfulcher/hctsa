function ord = BF_linkageOrdering(distMat,links)
% BF_linkageOrdering attempts to use optimalleaforder for dendrogram orderings

% ------------------------------------------------------------------------------
% Copyright (C) 2018, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

%-------------------------------------------------------------------------------
% Find the size threshold:
if any(size(distMat)==1)
    % A vector of pairwise distances
    numItems = (1+sqrt(1+8*length(distMat)))/2;
else
    % A square matrix of pairwise distances
    numItems = length(distMat);
end

%-------------------------------------------------------------------------------
% Do either optimalleaforder, or normal dendrogram ordering, applying a threhold
% on the number of items for which it is efficient to apply the optimalleaforder
% algorithm. 2000 is the number

if numItems < 2000 % small enough to try optimalleaforder
    try
        ord = optimalleaforder(links,distMat);
        % [~,~,ord] = dendrogram(links,0,'r',ord);
        fprintf(1,'Using optimalleaforder reordering!\n');
    catch
        fprintf(1,'Using dendrogram reordering.\n');
        [~,~,ord] = dendrogram(links,0);
    end
else
    fprintf(1,'Too many objects to reorder using optimalleaforder, using dendrogram instead.\n');
    [~,~,ord] = dendrogram(links,0);
end

end
