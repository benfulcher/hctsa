function [tsTime,masterTime,masterIDs] = TS_TotalCalcTime(whatData)
% TS_TotalCalcTime   Total calculation time per time series, counting each master
%                    operation once.
%
% TS_CalcTime stores, for every operation (feature), the time taken to evaluate
% its *master* operation -- so all operations drawing on the same master
% operation repeat the same time. Summing TS_CalcTime across operations therefore
% overcounts the true cost (by a factor of the number of outputs per master
% operation, ~8x on average for the full hctsa library). This function collapses
% TS_CalcTime to one time per master operation before summing.
%
%---INPUTS:
% whatData: an HCTSA .mat file name, or a structure loaded from one
%           (default: 'raw', i.e., HCTSA.mat)
%
%---OUTPUTS:
% tsTime: (numTimeSeries x 1) total calculation time (s) for each time series.
%         Master operations with no recorded time (e.g., errors, NaN outputs)
%         contribute nothing; a time series with no recorded times gets NaN.
% masterTime: (numTimeSeries x numMasters) calculation time (s) of each master
%             operation on each time series (NaN where not recorded).
% masterIDs: (numMasters x 1) the master operation ID of each column of masterTime.
%
%---EXAMPLE USAGE:
% tsTime = TS_TotalCalcTime('HCTSA.mat');
% fprintf(1,'Mean time per series: %s\n',BF_TheTime(mean(tsTime,'omitnan')));

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

if nargin < 1 || isempty(whatData)
    whatData = 'raw';
end
if ischar(whatData) && strcmp(whatData,'raw')
    whatData = 'HCTSA.mat';
end

TS_CalcTime = TS_GetFromData(whatData,'TS_CalcTime');
Operations = TS_GetFromData(whatData,'Operations');
if isempty(TS_CalcTime)
    error('No TS_CalcTime found in the data provided.');
end

% Group operations by master operation. Every operation of a master stores the
% same time (or NaN if it didn't return a good value), so take the max over each
% group, ignoring NaNs, to recover the master operation's time:
[masterIDs,~,groupInd] = unique(Operations.MasterID);
numMasters = length(masterIDs);
masterTime = nan(size(TS_CalcTime,1),numMasters);
for j = 1:numMasters
    masterTime(:,j) = max(TS_CalcTime(:,groupInd==j),[],2,'omitnan');
end

tsTime = sum(masterTime,2,'omitnan');
tsTime(all(isnan(masterTime),2)) = NaN;

end
