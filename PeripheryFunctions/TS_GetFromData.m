function result = TS_GetFromData(dataSource,dataField)
% TS_GetFromData   Load a given field from either a data file, or data structure
%                   (loaded from file).
%
%---INPUTS:
% dataSource: either a .mat filename or a structure generated from loading an
%               HCTSA file, e.g.,
%               dataSource = load('HCTSA.mat');
%               (the shorthands 'raw', 'loc', 'norm', and 'cl' are also accepted,
%               cf. TS_LoadData)
% dataField: the variable to extract
%
%---OUTPUTS:
% result: the variable extracted from either the file or structure (depending
%           on the type of dataSource). Empty if the variable is not present.
%           An error is thrown if the file does not exist.
%
% The placeholder returned by TS_LoadData for structure inputs,
% '--INPUT_STRUCTURE--', has no file to read from and always gives empty.
%
%---EXAMPLE USAGE:
% MasterOperations = TS_GetFromData('HCTSA_EEGDataset.mat','MasterOperations');

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

%-------------------------------------------------------------------------------
% In some cases, you provide a structure with the pre-loaded data already in it
% e.g., as a whatDataFile = load('HCTSA.mat');
%-------------------------------------------------------------------------------
if isstruct(dataSource)
    if isfield(dataSource,dataField)
        result = dataSource.(dataField);
    else
        result = [];
    end

%-------------------------------------------------------------------------------
% Often you provide a .mat file name: load in from file
%-------------------------------------------------------------------------------
elseif ischar(dataSource)

    % Decode shorthands for default filenames 'raw', 'loc', 'norm', and 'cl':
    switch dataSource
    case {'raw','loc'} % the raw, un-normalized data:
        dataSource = 'HCTSA.mat';
    case 'norm' % the normalized data:
        dataSource = 'HCTSA_N.mat';
    case 'cl' % the clustered data:
        dataSource = 'HCTSA_N.mat';
    end

    if strcmp(dataSource,'--INPUT_STRUCTURE--')
        % Name given by TS_LoadData to a pre-loaded structure: nothing to read
        result = [];
    else
        if ~exist(dataSource,'file')
            error('%s not found',dataSource);
        end
        % Only load the variable if it is in the file (avoids load's warning)
        fileVarsStruct = whos('-file',dataSource);
        if ismember(dataField,{fileVarsStruct.name})
            loadAgain = load(dataSource,dataField);
            result = loadAgain.(dataField);
        else
            result = [];
        end
    end
else
    error('Unknown data source type: %s',class(dataSource));
end

% Check for legacy structure array format and convert to table:
if ismember(dataField,{'TimeSeries','Operations','MasterOperations'}) && isstruct(result)
    warning('Loaded metadata, %s, is still in the old structure-array format; converted to table',dataField)
    result = struct2table(result);
end

end
