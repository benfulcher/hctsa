function OutputToCSV(whatData,writeTimeSeriesData)
% OutputToCSV     Outputs data to csv for analysis in other enrivonments
%
%---INPUTS:
% whatData, which HCTSA.mat file to use (default: HCTSA.mat)
% writeTimeSeriesData, (logical) whether to also output time-series data to file

% ------------------------------------------------------------------------------
% Copyright (C) 2017, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
% Check Inputs:
%-------------------------------------------------------------------------------
if nargin < 1
    whatData = 'raw';
end
if nargin < 2
    writeTimeSeriesData = false;
end
%-------------------------------------------------------------------------------

% Load in the data:
[TS_DataMat,TimeSeries,Operations,whatData] = TS_LoadData(whatData);

% Get the quality info:
TS_Quality = TS_GetFromData(whatData,'TS_Quality');

% Put NaNs in TS_DataMat where good quality is lacking (ignoring the code):
TS_DataMat(TS_Quality~=0) = NaN;

%-------------------------------------------------------------------------------
% Output data matrix to file:
%-------------------------------------------------------------------------------
fileName = 'hctsa_datamatrix.csv';
dlmwrite(fileName,TS_DataMat,'delimiter',',');
fprintf(1,'Wrote data to %s\n',fileName);

%-------------------------------------------------------------------------------
% Output time-series info to file:
%-------------------------------------------------------------------------------
combinedStrings = arrayfun(@(x)sprintf('%s (%s)',TimeSeries(x).Name,...
        TimeSeries(x).Keywords),1:length(TimeSeries),'UniformOutput',false)';
% Write it out
fileName = 'hctsa_timeseries-info.csv';
fid = fopen(fileName,'w');
for i = 1:length(TimeSeries)
    fprintf(fid,'%s,%s\n',TimeSeries(i).Name,TimeSeries(i).Keywords);
end
fclose(fid);
fprintf(1,'Wrote time-series info to %s\n',fileName);

%-------------------------------------------------------------------------------
% Output feature info to file:
%-------------------------------------------------------------------------------
fileName = 'hctsa_features.csv';
fid = fopen(fileName,'w');
for i = 1:length(Operations)
    fprintf(fid,'%s,%s\n',Operations(i).Name,Operations(i).CodeString);
end
fclose(fid);
fprintf(1,'Wrote feature info to %s\n',fileName);

%-------------------------------------------------------------------------------
% Output time-series data to file:
%-------------------------------------------------------------------------------
if writeTimeSeriesData
    fprintf(1,'Writing time-series data for %u time series...\n',length(TimeSeries));
    fileName = 'hctsa_timeseries-data.csv';
    fid = fopen(fileName,'w');
    for i = 1:length(TimeSeries)
        x = TimeSeries(i).Data;
        L = length(x);
        for t = 1:L-1
            fprintf(fid,'%.6g,',x(t)); % up to last value with commas
        end
        fprintf(fid,'%.6g\n',x(L)); % final value with new line
    end
    fclose(fid);
    fprintf(1,'Wrote time-series data to %s\n',fileName);
end

end
