function TS_ConvertToTable(hctsaFile)
% TS_ConvertToTable     Converts a file using the old structure array format -> table
%
%---INPUTS:
% hctsaFile, an hctsa file to convert

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

if nargin < 1
    hctsaFile = 'HCTSA.mat';
end

[~,TimeSeries,Operations,whatDataFile] = TS_LoadData(hctsaFile,false);
MasterOperations = TS_GetFromData(whatDataFile,'MasterOperations');

% Convert the three structure arrays to tables:
if isstruct(TimeSeries)
    TimeSeries = struct2table(TimeSeries);
end
if isstruct(Operations)
    Operations = struct2table(Operations);
end
if isstruct(MasterOperations)
    MasterOperations = struct2table(MasterOperations);
end

% Save back to file
fprintf(1,'Converted TimeSeries, Operations, and MasterOperations to table, saving back to %s...',whatDataFile);
save(whatDataFile,'TimeSeries','Operations','MasterOperations','-append')
fprintf(1,'Saved.\n');

end
