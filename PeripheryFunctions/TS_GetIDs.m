function IDs = TS_GetIDs(theFields,whatData,tsOrOps,nameOrKeywords)
% TS_GetIDs   Retrieve IDs of time series (or operations) in an hctsa dataset
%               using that matches a string in either the Name or Keyword
%               fields.
%
%---INPUTS:
% theFields, the fields to match on (string).
% whatData, the source of the hctsa dataset (e.g., a filename, cf. TS_LoadData).
%           (default: 'norm')
% tsOrOps, whether to retrieve IDs for TimeSeries ('ts', default) or
%           Operations ('ops').
% nameOrKeywords, what field to match, either the "Name" field or the
% "Keyword" field
%
%---OUTPUTS:
% IDs, a (sorted) vector of IDs matching the field constraint provided.
%
%---EXAMPLE USAGE:
% >> ts_IDs = TS_GetIDs('noisy','norm','ts');
% This retrieves the IDs of time series in 'HCTSA_N.mat' (specifying 'norm') that
% contain the keyword 'noisy'.
%
% >> op_IDs = TS_GetIDs('entropy','norm','ops');
% This retrieves the IDs of operations in HCTSA_N.mat that have been tagged
% with the keyword 'entropy'.

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
% Check inputs, set defaults:
%-------------------------------------------------------------------------------
if nargin < 2
    whatData = 'norm';
end
if nargin < 3
    tsOrOps = 'ts';
end
if nargin < 4
    nameOrKeywords = 'Keywords';
end

%-------------------------------------------------------------------------------
% Load data:
%-------------------------------------------------------------------------------
[~,TimeSeries,Operations,theDataFile] = TS_LoadData(whatData);

%-------------------------------------------------------------------------------
% Match time series/operations on an input keyword:
%-------------------------------------------------------------------------------
switch tsOrOps
case 'ts'
    theDataTable = TimeSeries;
case 'ops'
    theDataTable = Operations;
otherwise
    error('Specify ''ts'', ''ops'', or ''opsName''');
end

% OC: below has been changed so that we return a set in the same order as
%       above, and nan's otherwise
matches = nan(length(theFields),1);
for i = 1:length(theFields)
    cmatch = find(strcmp(theDataTable.(nameOrKeywords),theFields{i}));
    if ~isempty(cmatch)
        matches(i) = cmatch;
    end
end

foundIDs = ~isnan(matches);

IDs = zeros(size(matches));
IDs(foundIDs) = theDataTable.ID(matches(foundIDs));
IDs(~foundIDs) = nan;

if all(isnan(IDs))
    warning('No matches to ''%s'' found in %s',theFields,theDataFile)
end

end
