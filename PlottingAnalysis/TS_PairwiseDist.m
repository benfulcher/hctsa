% --------------------------------------------------------------------------
% TS_PairwiseDist
% --------------------------------------------------------------------------
% 
% Computes pairwise distances between all objects (either TimeSeries or
% Operations) and stores the data back in the local .mat file.
% 
% Warning: this can be large.
% 
%---INPUTS:
%
% tsOrOps: Compute pairwise distances between all pairs of time series 
%               ('ts', default), or operations ('ops')
%               
% whatData: 'orig' (load data from HCTSA_loc.mat)
%           'norm' (default: load data from HCTSA_N.mat)
%
% distanceMetric: what distance metric to use, e.g., 'euclidean' (default for
%                 time series) or 'corr_fast' (default for operations)
%
% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function Dij = TS_PairwiseDist(tsOrOps,whatData,distanceMetric)
    
% ------------------------------------------------------------------------------
% Check inputs:
% ------------------------------------------------------------------------------

if nargin < 1 || isempty(tsOrOps)
    tsOrOps = 'ts';
    fprintf(1,'Time series by default.\n');
end
if ~ismember(tsOrOps,{'ts','ops'})
    error('Invalid specifier ''%s'' (should be ''ts'' or ''ops'')',tsOrOps);
end

if nargin < 2
    whatData = 'norm';
end

if nargin < 3
    switch tsOrOps
    case 'ts'
        distanceMetric = 'euclidean';
    case 'ops'
        distanceMetric = 'corr_fast';
    end
end

% ------------------------------------------------------------------------------
% Load data
% ------------------------------------------------------------------------------
switch whatData
    case 'orig'
        theFile = 'HCTSA_loc.mat';
    case 'norm'
        theFile = 'HCTSA_N.mat';
    otherwise
        error('Invalid specifier ''%s''.',whatData);
end
load(theFile,'TS_DataMat');
   
% ------------------------------------------------------------------------------
% Compute pairwise distances between all objects
% ------------------------------------------------------------------------------ 
if strcmp(tsOrOps,'ops')
    TS_DataMat = TS_DataMat';
end

fprintf(1,'Computing %s distances between %u objects...',distanceMetric,size(TS_DataMat,1));
Dij = BF_pdist(TS_DataMat,distanceMetric,1);
fprintf(1,' Done.\n');

% ------------------------------------------------------------------------------
% Save back
% ------------------------------------------------------------------------------
fprintf(1,'Saving back to %s...',theFile);
switch tsOrOps
case 'ts'
    distanceMetricTimeSeries = distanceMetric;
    Dij_TimeSeries = Dij;
    save(theFile,'Dij_TimeSeries','distanceMetricTimeSeries','-append')
case 'ops'
    distanceMetricOperations = distanceMetric;
    Dij_Operations = Dij;
    save(theFile,'Dij_Operations','distanceMetricOperations','-append')
end
fprintf(1,' Done.\n');

end