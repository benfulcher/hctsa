function [IDs,notIDs] = TS_GetIDs(theMatchString,whatData,tsOrOps,nameOrKeywords)
% TS_GetIDs   Retrieve IDs of time series (or operations) in an hctsa dataset
%               using that matches a string in either the Name or Keyword
%               fields.
%
%---INPUTS:
% theMatchString, the string to match
% whatData, the source of the hctsa dataset (e.g., a filename, cf. TS_LoadData).
%           (default: 'norm')
% tsOrOps, whether to retrieve IDs for TimeSeries ('ts', default) or
%           Operations ('ops')
% nameOrKeywords, what field to match: 'Keywords' (default) or 'Name'
%
%---OUTPUTS:
% IDs, a (sorted) vector of IDs matching the field constraint provided.
%
%---EXAMPLE USAGE:
% >> tsIDs = TS_GetIDs('noisy','norm','ts','Keywords');
% This retrieves the IDs of time series in 'HCTSA_N.mat' (specifying 'norm') that
% contain the keyword 'noisy'.
%
% >> opIDs = TS_GetIDs('entropy','norm','ops','Keywords');
% This retrieves the IDs of operations in 'HCTSA_N.mat' that have been tagged
% with the keyword 'entropy'.
%
% NAME MODE:
% >> opIDs = TS_GetIDs('length','norm','ops','Name');
% >> opIDs = TS_GetIDs({'length','mean'},'norm','ops','Name');

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
if istable(whatData)
    theDataTable = whatData;
    theDataFile = 'the table provided';
else
    switch tsOrOps
    case 'ts'
        % Retrieve the TimeSeries table:
        [~,theDataTable,~,theDataFile] = TS_LoadData(whatData);
    case 'ops'
        % Retrieve the Operations table:
        [~,~,theDataTable,theDataFile] = TS_LoadData(whatData);
    otherwise
        error('Specify ''ts'' (time series) or ''ops'' (operations)');
    end
end

%-------------------------------------------------------------------------------
%------------------------------------------------------------------------------
switch nameOrKeywords
    %--------------------------------------------------------------------------
    % KEYWORDS: match IDs of time series or operations by exact keyword match
    %--------------------------------------------------------------------------
    case {'keywords','Keywords'}
        % The cell of comma-delimited keyword strings:
        theKeywordCell = theDataTable.Keywords;

        % Split into sub-cells using comma delimiter:
        Keywords = SUB_cell2cellcell(theKeywordCell);

        % Find objects with a keyword that matches the input string:
        matches = cellfun(@(x)any(ismember(theMatchString,x)),Keywords);

        % Return the IDs of the matches:
        IDs = theDataTable.ID(matches);

        % Check for empty:
        if isempty(IDs)
            warning('No matches to ''%s'' found in %s',theMatchString,theDataFile)
        end

        % Also provide IDs not matching the constraint, if required
        if nargout > 1
            notIDs = setxor(IDs,theDataTable.ID);
        end

    case {'name','Name'}
        %----------------------------------------------------------------------
        % NAME: Find match for each element of input, and return an ordered set
        % (with NaN when we don't find a match)
        %----------------------------------------------------------------------

        assert(nargout == 1,'One output argument (IDs) required for ''Name'' input.')

        if ischar(theMatchString)
            theMatchString = {theMatchString};
        end
        numStrings = length(theMatchString);
        matches = nan(numStrings,1);
        for i = 1:numStrings
            cmatch = find(strcmp(theDataTable.Name,theMatchString{i}));
            if ~isempty(cmatch)
                matches(i) = cmatch;
            end
        end

        foundIDs = ~isnan(matches);

        IDs = zeros(size(matches));
        IDs(foundIDs) = theDataTable.ID(matches(foundIDs));
        IDs(~foundIDs) = nan;

        if all(isnan(IDs))
            warning('No matches to the %u input strings found in %s',numStrings,theDataFile)
        end

    otherwise
        error('Must specify either ''Keywords'' or ''Name''');
end

end
