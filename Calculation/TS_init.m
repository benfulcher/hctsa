function TS_init(INP_ts,INP_mops,INP_ops,beVocal,outputFile)
% TS_init       Takes in time series, master operation, and operation input
% files and produces a formatted HCTSA .mat file
%
% This function is used instead to run hctsa analysis without a linked mySQL database.
%
%---EXAMPLE USAGE:
% Initiate an HCTSA.mat file on a custom time-series dataset using default
% feature library, using a formatted input file, 'my_TS_INP_file.mat'
% >> TS_init('my_TS_INP_file.mat');
%
%---INPUTS:
% INP_ts: A time-series input file
% INP_mops: A master operations input file
% INP_ops: An operations input file
% beVocal: Whether to display details of the progress of the script to screen.
%           a 3-vector, specifying for 1. time series, 2. master operations,
%           and 3. operations.
% outputFile: Specify an alternative output filename
%
%---OUTPUTS:
% Writes output into HCTSA.mat (or specified custom filename)

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

% ------------------------------------------------------------------------------
% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 1 || isempty(INP_ts)
    error('Please supply a formatted time-series input file (see documentation for details).');
end
if nargin < 2 || isempty(INP_mops)
    INP_mops = 'INP_mops.txt';
end
if nargin < 3 || isempty(INP_ops)
    INP_ops = 'INP_ops.txt';
end
if nargin < 4
    if nargin < 2
        beVocal = [1,0,0]; % by default helps you just for the time series input file you provided
    elseif nargin < 3
        beVocal = [1,1,0]; % by help you through the master operations too
    else
        beVocal = [1,1,1]; % Provided all custom input files--walks you through all of them
    end
end
if length(beVocal)==1
    beVocal = ones(3,1)*beVocal;
end
if nargin < 5
    outputFile = 'HCTSA.mat';
end

% ------------------------------------------------------------------------------
% First check if you're about to overwrite an existing file
% ------------------------------------------------------------------------------
if exist(['./',outputFile],'file')
    reply = input(sprintf(['Warning: %s already exists -- if you continue, this ' ...
        'file will be overwritten.\n[press ''y'' to continue]'],outputFile),'s');
    if ~strcmp(reply,'y')
        return
    end
end

% ------------------------------------------------------------------------------
% Get time series, operations, master operations into structure arrays
% ------------------------------------------------------------------------------
isEmptyStruct = @(x) length(fieldnames(x))==0;

TimeSeries = SQL_add('ts', INP_ts, 0, beVocal(1));
if isEmptyStruct(TimeSeries), return; end % The user did not approve of the set of inputs
numTS = length(TimeSeries);

MasterOperations = SQL_add('mops', INP_mops, 0, beVocal(2))';
if isEmptyStruct(MasterOperations), return; end % The user did not approve of the set of inputs
numMops = length(MasterOperations);

Operations = SQL_add('ops', INP_ops, 0, beVocal(3));
if isEmptyStruct(Operations), return; end % The user did not approve of the set of inputs
numOps = length(Operations);

%-------------------------------------------------------------------------------
% Link operations to their masters using label matching
% and update the structure arrays using the TS_LinkOperationsWithMasters function
%-------------------------------------------------------------------------------

[Operations, MasterOperations] = TS_LinkOperationsWithMasters(Operations,MasterOperations);

% MasterOperations may have been trimmed by TS_LinkOperationsWithMasters:
numMops = length(MasterOperations);

% ------------------------------------------------------------------------------
% Generate the TS_DataMat, TS_Quality, and TS_CalcTime matrices
% ------------------------------------------------------------------------------
% All NaNs -> NULL (haven't yet been calculated)
TS_DataMat = nan(numTS,numOps);
TS_Quality = nan(numTS,numOps);
TS_CalcTime = nan(numTS,numOps);

%-------------------------------------------------------------------------------
% Get git information to keep track of the version of code used at the time of TS_init
%-------------------------------------------------------------------------------
gitInfo = TS_AddGitInfo();

% ------------------------------------------------------------------------------
% Save to file
% ------------------------------------------------------------------------------
% Set a flag, fromDatabase, that tells you that you that this was generated by
% TS_init and shouldn't be written back to a database
fromDatabase = 0;
save(outputFile,'TimeSeries','Operations','MasterOperations',...
            'TS_DataMat','TS_Quality','TS_CalcTime','fromDatabase','gitInfo','-v7.3');

fprintf(1,'Successfully initialized %s with %u time series, %u master operations, and %u operations\n',...
                        outputFile,numTS,numMops,numOps);

end
