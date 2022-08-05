function [TS_DataMat,TimeSeries] = FilterLabeledTimeSeries(TS_DataMat,TimeSeries)
%
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

% Check that time series have been labeled
if ~ismember('Group',TimeSeries.Properties.VariableNames)
    error('Time series are not yet labeled. Please run TS_LabelGroups on your hctsa dataset.');
end

% Check that labeling is done with categoricals (rather than legacy integer labels)
if isnumeric(TimeSeries.Group)
    error(['Time-series group labels are assigned as numeric labels (old format).',...
                    ' Please re-run TS_LabelGroups on your hctsa dataset.']);
end
isLabeled = ~isundefined(TimeSeries.Group);
if any(~isLabeled)
    TimeSeries = TimeSeries(isLabeled,:);
    TS_DataMat = TS_DataMat(isLabeled,:);
    fprintf(1,'Only considering the %u time series that have labels.\n',sum(isLabeled));
end

end
