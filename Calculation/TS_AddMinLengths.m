function [Operations,MasterOperations] = TS_AddMinLengths(Operations,MasterOperations,minLengthFile)
% TS_AddMinLengths  Attach minimum time-series lengths to a feature set
%
% Adds a MinLength column to the MasterOperations and Operations tables. A
% master operation is not evaluated on a time series shorter than its
% MinLength: all of its operations are then NaN (TS_Quality = 2, 'not
% appropriate'), the outcome hctsa uses for an algorithm that cannot be
% applied to data like this. An operation's own MinLength sets that single
% output to NaN on shorter series while the master operation's other outputs
% are still computed. A MinLength of 0 means no bound.
%
% Bounds refer to the length of the series the master operation is applied
% to, which is the decimated series for master operations that use
% BF_PreProcess(x_z,'decimate_ac1e') (see TS_CalculateFeatureVector).
%
%---INPUTS:
% Operations, MasterOperations: tables of operations and master operations,
%           linked by TS_LinkOperationsWithMasters (Operations.MasterID).
% minLengthFile: a minimum-length input file (e.g., 'INP_minlength_hctsa.txt'),
%           or a feature-set name (e.g., 'hctsa'), or empty for no bounds.
%           Each line of the file holds a key and a minimum length, separated by
%           whitespace; lines starting with # are comments. The key is either a
%           master-operation label (e.g., SP_Summaries_welch_rect), giving the
%           bound for the master operation, or an operation code string
%           (MasterLabel.field, e.g., ST_RangeEvolve.p1), giving a bound for that
%           operation alone. Keys that match nothing in the tables are ignored, so
%           one file serves any subset of its feature set.
%
%---OUTPUTS:
% Operations, MasterOperations: the input tables, each with a MinLength column.

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

if nargin < 3
    minLengthFile = '';
end
if ~ismember('MasterID',Operations.Properties.VariableNames)
    error('Link Operations to MasterOperations (TS_LinkOperationsWithMasters) before adding minimum lengths.');
end

MasterOperations.MinLength = zeros(height(MasterOperations),1);
Operations.MinLength = zeros(height(Operations),1);

% A feature-set name stands for that feature set's minimum-length file:
if ~isempty(minLengthFile) && ~endsWith(minLengthFile,'.txt')
    [~,~,minLengthFile] = TS_GiveMeInputFiles(minLengthFile);
end
if isempty(minLengthFile)
    return
end
if ~exist(minLengthFile,'file')
    error('Unknown minimum-length file: ''%s''',minLengthFile);
end

fid = fopen(minLengthFile);
if fid == -1
    error('Could not open the minimum-length file ''%s''',minLengthFile);
end
dataIn = textscan(fid,'%s%f','CommentStyle','#');
fclose(fid);
keys = dataIn{1};
minLengths = dataIn{2};
if numel(minLengths) ~= numel(keys) || any(isnan(minLengths)) || ...
        any(minLengths < 0 | minLengths ~= round(minLengths))
    error('Each line of %s must be a key followed by a non-negative integer.',minLengthFile);
end
if numel(unique(keys)) < numel(keys)
    error('%s lists a key more than once.',minLengthFile);
end

% Master operations (keys without a '.'):
isOpKey = contains(keys,'.');
[isMatch,ind] = ismember(MasterOperations.Label,keys(~isOpKey));
masterBounds = minLengths(~isOpKey);
MasterOperations.MinLength(isMatch) = masterBounds(ind(isMatch));

% Single operations (keys of the form MasterLabel.field), matched on their code string:
[isMatch,ind] = ismember(Operations.CodeString,keys(isOpKey));
opBounds = minLengths(isOpKey);
Operations.MinLength(isMatch) = opBounds(ind(isMatch));

end
