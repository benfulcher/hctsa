function [distVec,links,distVec0,numItems] = BF_DoLinkage(distMat,whatDistance,linkageMeth)
% BF_DoLinkage     Process an input distance matrix for linkage clustering
% (including custom processing for absolute correlation distances)

%---INPUTS:
%
%
%---OUTPUTS:


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
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%% Check inputs
%-------------------------------------------------------------------------------
if nargin < 2
    whatDistance = 'corr';
end
if nargin < 3
    linkageMeth = 'average';
end

%-------------------------------------------------------------------------------
%% Work with vector form of distances
% (hopefully enough memory)
%-------------------------------------------------------------------------------
if any(size(distMat)==1)
    distVec = distMat;
else
    distVec = squareform(distMat); % convert for vector version
end
% Solve for number of items from quadratic formula:
numItems = (1 + sqrt(1 + 8*length(distVec)))/2;

%-------------------------------------------------------------------------------
%% Convert to absolute correlations:
%-------------------------------------------------------------------------------
distVec0 = distVec; % keep distVec0 as the input distance matrix in all cases
% distVec changes only for abscorr options
if ismember(whatDistance,{'abscorr','abscorr_ii'})
    % Compute distances on absolute correlation distances (where sign of correlation
    % is irrelevant, it's the magnitude that's important):
    distVec = 1 - abs(1 - distVec);
end

% We need the matrix form also:
distMat = squareform(distVec);

%-------------------------------------------------------------------------------
%% Do the linkage clustering:
%-------------------------------------------------------------------------------
fprintf(1,'Computing linkage information for %u x %u data using %s clustering...',...
            numItems,numItems,linkageMeth);
links = linkage(distVec,linkageMeth);
fprintf(1,' Done.\n');

end
