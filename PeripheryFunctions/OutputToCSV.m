function OutputToCSV(whatData,writeTimeSeriesData,writeMasterFeatures)
% OutputToCSV     Outputs data to csv for external analysis
%
%---INPUTS:
% whatData, which hctsa .mat file to use (default: HCTSA.mat)
% writeTimeSeriesData, (logical) whether to also output time-series data to file
% writeMasterFeatures, (logical) whether to also output master operation info

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
% Check inputs and set defaults
%-------------------------------------------------------------------------------
if nargin < 1 || isempty(whatData)
    whatData = 'raw';
end
if nargin < 2 || isempty(writeTimeSeriesData)
    writeTimeSeriesData = false;
end
if nargin < 3
    writeMasterFeatures = false;
end
theDelimiter = ',';
%-------------------------------------------------------------------------------

% Load in the hctsa data:
[TS_DataMat,TimeSeries,Operations,whatData] = TS_LoadData(whatData);

% Get the master operations
if writeMasterFeatures
    MasterOperations = TS_GetFromData(whatData,'MasterOperations');
end

% Get the quality info:
TS_Quality = TS_GetFromData(whatData,'TS_Quality');

% Put NaNs in TS_DataMat where good quality is lacking (ignoring the code):
TS_DataMat(TS_Quality~=0) = NaN;

%-------------------------------------------------------------------------------
% Output data matrix to file:
%-------------------------------------------------------------------------------
fileName = 'hctsa_datamatrix.csv';
dlmwrite(fileName,TS_DataMat,'delimiter',theDelimiter);
fprintf(1,'Wrote feature matrix to %s.\n',fileName);

%-------------------------------------------------------------------------------
% Output time-series info to .csv file:
%-------------------------------------------------------------------------------
fileName = 'hctsa_timeseries-info.csv';
if strcmp(theDelimiter,',')
    doQuoteStrings = true;
    % warning(['The comma-delimited Keywords variable may cause trouble for ',...
        % 'the comma-delimited %s file: quoting strings should help'],fileName);
else
    doQuoteStrings = false;
end

% We don't want the data to be written out in this info file:
v = version;
if str2double(v(1:3)) < 9.4 
    % 2017 or older MATLAB version
    TimeSeries(:,'Data') = [];
else
    TimeSeries = removevars(TimeSeries,{'Data'});
end
writetable(TimeSeries,fileName,'FileType','text','WriteVariableNames',true,...
                'Delimiter',theDelimiter,'QuoteStrings',doQuoteStrings)
fprintf(1,'Wrote the TimeSeries table to %s.\n',fileName);

%-------------------------------------------------------------------------------
% Output feature info to file:
%-------------------------------------------------------------------------------
fileName = 'hctsa_features.csv';
writetable(Operations,fileName,'FileType','text','WriteVariableNames',true,...
                            'Delimiter',theDelimiter)
fprintf(1,'Wrote the Operations table to %s.\n',fileName);

%-------------------------------------------------------------------------------
% Output time-series data to file:
%-------------------------------------------------------------------------------
if writeTimeSeriesData
    numTimeSeries = height(TimeSeries);
    fprintf(1,'Writing time-series data for %u time series...\n',numTimeSeries);
    fileName = 'hctsa_timeseries-data.csv';
    fid = fopen(fileName,'w');
    for i = 1:numTimeSeries
        x = TimeSeries.Data{i};
        L = length(x);
        for t = 1:L-1
            fprintf(fid,'%.6g,',x(t)); % up to last value with commas
        end
        fprintf(fid,'%.6g\n',x(L)); % final value with new line
    end
    fclose(fid);
    fprintf(1,'Wrote comma-delimited time-series data (6-decimal precision) to %s.\n',fileName);
end

%-------------------------------------------------------------------------------
% Output master features info to file:
%-------------------------------------------------------------------------------
if writeMasterFeatures
    fileName = 'hctsa_masterfeatures.csv';
    writetable(MasterOperations,fileName,'FileType','text','WriteVariableNames',true,...
                                'Delimiter',theDelimiter)
    fprintf(1,'Wrote the MasterOperations table to %s.\n',fileName);
end

%-------------------------------------------------------------------------------
% Output time-series group info to file:
% (it's accessible from the time-series table)
%-------------------------------------------------------------------------------
% if ismember('Group',TimeSeries.Properties.VariableNames)
%     fileName = 'hctsa_grouplabel-info.csv';
%     fid = fopen(fileName,'w');
%     % Header:
%     fprintf(fid,'%s%s%s\n','Name',theDelimiter,'Group');
%     for i = 1:height(TimeSeries)
%         fprintf(fid,'%s%s%u\n',TimeSeries.Name{i},theDelimiter,TimeSeries.Group(i));
%     end
%     fclose(fid);
%     fprintf(1,'Wrote time-series info to %s.\n',fileName);
% end


end
