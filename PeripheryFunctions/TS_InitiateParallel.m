function success = TS_InitiateParallel(doInitiate)
% Initiates a parallel session for hctsa
%
%---INPUT:
% doInitiate: = 1 for the first time you want to initialize all workers
%
%---OUTPUT:
% success -- = 1 if a matlabpool was successfully initiated (or if one already
%               exists)

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

% Check inputs:
if nargin < 1
    doInitiate = false;
end

% Check a license is available:
try
    BF_CheckToolbox('distrib_computing_toolbox')
catch
    fprintf(1,['License for Parallel Computing Toolbox could not be initiated' ...
                    ' -- cannot perform computations across multiple cores\n']);
    success = false;
    return
end
success = true;

%-------------------------------------------------------------------------------
matlabVersion = version('-release');

% Syntax changed in Matlab 2015a
try
    if str2num(matlabVersion(1:4)) >= 2014

        poolObj = gcp('nocreate'); % If no pool already, create a new one

        if isempty(poolObj) % no matlab pool started yet
            % Open pool of workers:
            poolObj = parpool;
            % Get number of workers:
            numWorkers = poolObj.NumWorkers;
            % User feedback:
            fprintf(1,['Matlab parallel processing pool opened with %u ' ...
                                    'workers\n'],numWorkers);
            % Regardless of what you input -- must initiate in this case:
            doInitiate = true;
        else
            % Get number of workers:
            numWorkers = poolObj.NumWorkers;
            fprintf(1,['Matlab parallel processing pool already open with ' ...
                                        '%u workers\n'],numWorkers);
        end
    else
        if (matlabpool('size') == 0)
            % Open pool of workers:
            matlabpool open;
            fprintf(1,['Matlab parallel processing pool opened with %u ' ...
                                    'workers\n'],matlabpool('size'));
        else
            fprintf(1,['Matlab parallel processing pool already open with ' ...
                                        '%u workers\n'],matlabpool('size'));
        end
    end
catch emsg
    warning('\nError starting parallel processing pool -- running serially instead:\n%s',emsg.message)
    success = false;
end

%-------------------------------------------------------------------------------
% Initiate hctsa stuff on all workers using pctRunOnAll
%-------------------------------------------------------------------------------
if doInitiate
    pctRunOnAll AddJavaClassPathInfoDynamics()
end

end
