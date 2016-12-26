function [Operations, MasterOperations] = TS_LinkOperationsWithMasters(Operations,MasterOperations)
% TS_LinkOperationsWithMasters    Link Operations with MasterOperations using Label field

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

numOps = length(Operations);
numMops = length(MasterOperations);

% ------------------------------------------------------------------------------
% Match operations to a master ID
% ------------------------------------------------------------------------------
for i = 1:numOps
    theMasterMatch = strcmp(Operations(i).Label,{MasterOperations.Label});
    if sum(theMasterMatch)==0
        error('No master match for operation: %s',Operations(i).Name);
    end
    Operations(i).MasterID = MasterOperations(theMasterMatch).ID;
end

% No longer need the label field
Operations = removeField(Operations,'Label');

%-------------------------------------------------------------------------------
% Check that all master operations are required
%-------------------------------------------------------------------------------
mastersNeeded = ismember([MasterOperations.ID],[Operations.MasterID]);
if ~all(mastersNeeded)
    warning(['%u/%u master operations are not used by the %u operations' ...
                         ' and will be removed.'],...
                        sum(~mastersNeeded),numMops,numOps);
    MasterOperations = MasterOperations(mastersNeeded);
end


% ------------------------------------------------------------------------------
function newStructArray = removeField(oldStructArray,fieldToRemove)

    theFieldNames = fieldnames(oldStructArray);

    fieldInd = strcmp(theFieldNames,fieldToRemove);

    oldCell = squeeze(struct2cell(oldStructArray));

    newCell = oldCell(~fieldInd,:);

    newStructArray = cell2struct(newCell,theFieldNames(~fieldInd));
end
% ------------------------------------------------------------------------------

end
