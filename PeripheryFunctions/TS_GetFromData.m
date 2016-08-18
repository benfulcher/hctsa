function result = TS_GetFromData(dataSource,dataField)
% TS_GetFromData   Load a given field from either a data file, or data structure
%                   (loaded from file).
%
%---INPUTS:
% dataSource: either a .mat filename or a structure generated from e.g.,
%               dataSource = load('HCTSA.mat');
% dataField: the variable to extract
%
%---OUTPUTS:
% Result is the variable extracted from either the file or structure (depending
% on the type of dataSource)

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
end

%-------------------------------------------------------------------------------
% Often you provide a .mat file name: load in from file
%-------------------------------------------------------------------------------
if ischar(dataSource)
    fileVarsStruct = whos('-file',dataSource);
    fileVars = {fileVarsStruct.name};
    if ismember('groupNames',fileVars)
        loadAgain = load(dataSource,dataField);
        result = loadAgain.(dataField);
    else
        result = [];
    end
end

end
