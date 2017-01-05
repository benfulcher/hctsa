function OutputToCSV(whatData)
% OutputToCSV     Outputs data to csv for analysis in other enrivonments
%
%---INPUTS:
% whatData, which HCTSA.mat file to use (default: HCTSA.mat)

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

if nargin < 1
    whatData = 'raw';
end

% Load in the data:
[TS_DataMat,TimeSeries,Operations,whatData] = TS_LoadData(whatData);

% Get the quality info:
load(whatData,'TS_Quality');

% Put NaNs in TS_DataMat where good quality is lacking (ignoring the code):
TS_DataMat(TS_Quality~=0) = NaN;

% Output data matrix to file:
dlmwrite('hctsa_datamatrix.csv',TS_DataMat,'delimiter',',');

% Output time series info to file:
combinedStrings = arrayfun(@(x)sprintf('%s (%s)',TimeSeries(x).Name,...
        TimeSeries(x).Keywords),1:length(TimeSeries),'UniformOutput',false)';
% Write it out
fid = fopen('hctsa_timeseries.csv','w');
for i = 1:length(TimeSeries)
    fprintf(fid,'%s,%s\n',TimeSeries(i).Name,TimeSeries(i).Keywords);
end
fclose(fid);

% Output feature info to file:
fid = fopen('hctsa_features.csv','w');
for i = 1:length(Operations)
    fprintf(fid,'%s,%s\n',Operations(i).Name,Operations(i).CodeString);
end
fclose(fid);

end
