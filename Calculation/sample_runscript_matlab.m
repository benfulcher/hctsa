function sample_runscript_matlab(doParallelize,saveIncrement,fileName)
% SAMPLE_RUNSCRIPT_MATLAB  Template function for looping over a hctsa analysis when
%                           using local Matlab files (initializing with TS_init).
%
%---INPUTS:
% doParallelize, (binary) whether to parallelize computations over available
%                   CPUs using Matlab's Parallel Computing Toolbox (default: 1).
% saveIncrement, saves the results back to the hctsa file after computing for
%               this many time series (default: 5)
% fileName, the source of the HCTSA dataset (a .mat file)
%
%---EXAMPLE USAGE:
% sample_runscript_matlab(1,'HCTSA.mat');

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
% Set parameters for computation:
%-------------------------------------------------------------------------------
if nargin < 1
    % Set to 1 to
    doParallelize = 1;
end

if nargin < 2
    saveIncrement = 5;
end

if nargin < 3
    fileName = 'HCTSA.mat';
end

% ------------------------------------------------------------------------------
%% Preparations:
% ------------------------------------------------------------------------------
% Get IDs of all time series in the hctsa datafile:
[~,TimeSeries,Operations] = TS_LoadData(fileName);
tsIDs = [TimeSeries.ID];

% Set the IDs to compute at each iteration:
ID_inc = 1:saveIncrement:length(tsIDs)+1;

% Provide a quick message:
fprintf(1,['About to calculate across %u time series and %u operations ' ...
                'saving back to file every %u time series\n'], ...
                length(TimeSeries),length(Operations),saveIncrement);

%-------------------------------------------------------------------------------
% Running calculations:
%-------------------------------------------------------------------------------
% Loop across time series:
for i = 1:length(ID_inc)-1
	fprintf(1,'\n\n\nWe''re looking at Time series with IDs from %u--%u\n\n\n', ...
                            	tsIDs(ID_inc(i)),tsIDs(ID_inc(i+1)-1));

    % Compute any missing values for this range of time series, then save back:
    TS_compute(doParallelize,tsIDs(ID_inc(i):ID_inc(i+1)-1),[],'missing',fileName);
end
