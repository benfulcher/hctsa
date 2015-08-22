function TS_local_clear_remove(tsOrOps,idRange,doRemove,whatData)
% TS_local_clear_remove     Clears or removed data from local files
%
%---INPUTS:
% tsOrOps -- either 'ts' or 'ops' for whether to work with either time series or operations
% idRange -- a vector of the ts_ids or op_ids in the database to remove
% doRemove -- whether to remove entries (specify 1), or just clear their data (specify 0)
% whatData -- the data to load (cf. TS_LoadData)

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>
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

%-------------------------------------------------------------------------------
%% Preliminaries and input checking
%-------------------------------------------------------------------------------

if nargin < 1
	tsOrOps = 'ts';
end
switch tsOrOps
case 'ts'
    theWhat = 'TimeSeries';
case 'ops'
    theWhat = 'Operations';
otherwise
    error('Specify either ''ts'' or ''ops''')
end

% Must specify a set of time series
if nargin < 2 || min(size(idRange)) ~= 1
	error('Specify a range of IDs');
end

if nargin < 3 % doRemove
    error('You must specify whether to remove the %s or just clear their data results',theWhat)
end

if nargin < 4
    whatData = 'loc'; % normally want to clear data from the local store
end

% ------------------------------------------------------------------------------
%% Load data
% ------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,whatDataFile] = TS_LoadData(whatData);

%-------------------------------------------------------------------------------
% Match IDs to indices
%-------------------------------------------------------------------------------
switch tsOrOps
case 'ts'
    dataStruct = TimeSeries;
case 'ops'
    dataStruct = Operations;
otherwise
    error('Specify either ''ts'' or ''ops''');
end
IDs = [dataStruct.ID];
doThese = ismember(IDs,idRange);

if ~any(doThese)
    error('No matches to %s',tsOrOps);
end

% ------------------------------------------------------------------------------
%% Provide some user feedback
% ------------------------------------------------------------------------------
if (doRemove == 0) % clear data
    doWhat = 'Clear data from';
else
    doWhat = '*PERMENANTLY REMOVE*';
end

iThese = find(doThese);
for i = 1:sum(doThese)
    fprintf(1,'%s [%u] %s\n',doWhat,IDs(iThese(i)),dataStruct(iThese(i)).Name)
end

if (doRemove == 0) % clear data
    reply = input(sprintf(['**Preparing to clear all calculated data for %u %s.\n' ...
                                '[press any key to continue]'], ...
                                    sum(doThese),theWhat),'s');
elseif doRemove == 1
    reply = input(sprintf(['Preparing to REMOVE %u %s -- DRASTIC STUFF! ' ...
                                'I HOPE THIS IS OK?!\n[press any key to continue]'], ...
                                sum(doThese),theWhat),'s');
else
    error('Specify either (0 to clear), or (1 to remove)')
end

% ------------------------------------------------------------------------------
%% Check what to clear/remove
% ------------------------------------------------------------------------------
if doRemove
    % Need to actually remove metadata entries:
    if strcmp(tsOrOps,'ts')
        TimeSeries(doThese) = [];
    else
        Operations(doThese) = [];
    end
end

TS_DataMat = f_clear_remove(TS_DataMat,doThese,tsOrOps,doRemove);

save(whatDataFile,'TS_DataMat','TimeSeries','Operations','-append');

% Repeat for any other matrix data files:
varNames = whos('-file',whatDataFile);
varNames = {varNames.name};
if ismember('TS_Quality',varNames)
    load(whatDataFile,'TS_Quality')
    TS_Quality = f_clear_remove(TS_Quality,doThese,tsOrOps,doRemove);
    save(whatDataFile,'TS_Quality','-append');
end
if ismember('TS_CalcTime',varNames)
    load(whatDataFile,'TS_CalcTime')
    TS_CalcTime = f_clear_remove(TS_CalcTime,doThese,tsOrOps,doRemove);
    save(whatDataFile,'TS_CalcTime','-append');
end

fprintf(1,'Saved back to %s\n',whatDataFile)

%-------------------------------------------------------------------------------
function A = f_clear_remove(A,ind,tsOrOps,doRemove)
    if doRemove
        if strcmp(tsOrOps,'ts')
            A(ind,:) = [];
        else
            A(:,ind) = [];
        end
    else
        if strcmp(tsOrOps,'ts')
            A(ind,:) = NaN;
        else
            A(:,ind) = NaN;
        end
    end
end
%-------------------------------------------------------------------------------

end
