function TS_CheckTimeSeries(whatData,whatID,doParallel)
% TS_CheckTimeSeries checks calculation of feature values for a given time series
% Specify the hctsa data file, and the ID of the time series to check against
% the operations.
%
%---INPUTS:
% whatData, the data to use (cf. TS_LoadData)
% whatID, the ID of the time series to check

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

%-------------------------------------------------------------------------------
% Check Inputs
%-------------------------------------------------------------------------------
if nargin < 1
    whatData = 'HCTSA.mat';
end
if nargin < 2
    whatID = [];
end
if nargin < 3
    doParallel = 1;
end

%-------------------------------------------------------------------------------
% Load in data, and filter
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations] = TS_LoadData(whatData);
MasterOperations = TS_GetFromData(whatData,'MasterOperations');
TS_Quality = TS_GetFromData(whatData,'TS_Quality');
if isempty(MasterOperations) || isempty(TS_Quality)
    error('MasterOperations, TS_Quality not found in the data source');
end
if length(TimeSeries) == 0
    error('No time series?!');
end
if isempty(whatID)
    idMatch = 1; % just take the first one
else
    idMatch = [TimeSeries.ID]==whatID;
end
if ~any(idMatch)
    error('ts_id %u not found in data source',whatID);
end
TimeSeries = TimeSeries(idMatch);
TS_DataMat = TS_DataMat(idMatch,:)';
TS_Quality = TS_Quality(idMatch,:)';

f = figure('color','w');
plot([TimeSeries.Data],'k')
title(TimeSeries.Name,'interpreter','none')

%-------------------------------------------------------------------------------
% Run the check:
%-------------------------------------------------------------------------------
[featureVector,~,calcQuality] = TS_CalculateFeatureVector(TimeSeries,...
                                    doParallel,Operations,MasterOperations,1,1);

%-------------------------------------------------------------------------------
% Compare output
%-------------------------------------------------------------------------------
% Quality of outputs:
misMatch = (TS_Quality~=calcQuality);
if any(misMatch)
    fprintf(1,'%u mismatches in quality... :(\n',sum(misMatch));
    % ------------
    % Text output:
    fmisMatch = find(misMatch);
    misMatchOps = Operations(misMatch);
    for i = 1:length(misMatchOps)
        fprintf('[%u] %s (%s, %s): %u (file) -> %u (now)\n',misMatchOps(i).ID,...
            misMatchOps(i).CodeString,...
            MasterOperations([MasterOperations.ID]==misMatchOps(i).MasterID).Code,...
            misMatchOps(i).Keywords,...
            TS_Quality(fmisMatch(i)),calcQuality(fmisMatch(i)));
    end
else
    fprintf(1,'All quality codes still match...! :)\n',sum(misMatch));
end

% first real outputs:
didWork = (calcQuality == 0);
matchMargin = (featureVector - TS_DataMat);
matchMarginProp = (featureVector - TS_DataMat)./TS_DataMat;
noMatch = ((abs(matchMarginProp) > 0.1) & didWork); % more than 0.1% different

% ------------
% Text output:
% ------------
% Find those that are stochastic and wouldn't be expected to reproduce:
dataStruct = struct();
dataStruct.Operations = Operations; dataStruct.TimeSeries = []; dataStruct.TS_DataMat = [];
stochasticIDs = TS_getIDs('stochastic',dataStruct,'ops');
isStochastic = ismember([Operations.ID],stochasticIDs)';

for j = 1:2
    if j==1
        % Labeled as stochastic first:
        fprintf(1,'\n\nOPERATIONS LABELED AS ''STOCHASTIC'' would not be expected to reproduce:\n\n\n');
        indx = find(noMatch & isStochastic);
        [~,ix] = sort(abs(matchMarginProp(noMatch & isStochastic)),'descend');
    else
        fprintf(1,'\n\nNOW OPERATIONS NOT LABELED AS ''STOCHASTIC'' would expected reproducibility:\n\n');
        indx = find(noMatch & ~isStochastic);
        [~,ix] = sort(abs(matchMarginProp(noMatch & ~isStochastic)),'descend');
    end
    indx = indx(ix);
    for i = 1:length(indx)
        ind = indx(i);
        fprintf('[%u] %s (%s | %s) [diff=%g, %.1f%%] %.8g (file) -> %.8g (now).\n',...
                Operations(ind).ID,...
                Operations(ind).CodeString,...
                MasterOperations([MasterOperations.ID]==Operations(ind).MasterID).Code,...
                Operations(ind).Keywords,...
                featureVector(ind)-TS_DataMat(ind),...
                100*(featureVector(ind)-TS_DataMat(ind))/TS_DataMat(ind),...
                TS_DataMat(ind),...
                featureVector(ind));
    end
end

end
