% ------------------------------------------------------------------------------
% compile_tisean
% ------------------------------------------------------------------------------
% This script attempts to compile the TISEAN package, putting the binaries into
% the Toolboxes/TiseanBinaries directory
% 
%---HISTORY:
% Ben Fulcher, 2015-03-06
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------


% ------------------------------------------------------------------------------
% Check we're in the correct folder
% ------------------------------------------------------------------------------

currentDir = pwd;

% Path split using platform-dependent separator
if isunix
    weHere = regexp(currentDir,'/','split');
else
    weHere = regexp(currentDir,'\','split');
end

if ~strcmp(weHere{end},'Toolboxes')
    error('This code must be run in the ''Toolboxes'' directory of the HCTSA package...')
end

% Sweet. Toolbox path is:
ToolDir = [currentDir '/'];

% ------------------------------------------------------------------------------
% Run system commands to compile Tisean
% ------------------------------------------------------------------------------
% Move to TISEAN directory
cd('Tisean_3.0.1');
% Configure:

