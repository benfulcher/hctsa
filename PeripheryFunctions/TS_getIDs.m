function ts_ids = TS_getIDs(theKeyword,whatData)
% TS_getIDs   Retrieve IDs of time series based on keyword matching.

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
% Check inputs, set defaults:
%-------------------------------------------------------------------------------
if nargin < 2
    whatData = 'norm';
end

%-------------------------------------------------------------------------------
% Load data:
%-------------------------------------------------------------------------------
[~,TimeSeries] = TS_LoadData(whatData);

%-------------------------------------------------------------------------------
% Match time series on an input keyword:
%-------------------------------------------------------------------------------
Keywords = SUB_cell2cellcell({TimeSeries.Keywords}); % Split into sub-cells using comma delimiter
matches = cellfun(@(x)any(ismember(theKeyword,x)),Keywords);
ts_ids = [TimeSeries(matches).ID];

end
