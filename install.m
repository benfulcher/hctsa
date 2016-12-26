% INSTALL   Installs the hctsa code package from scratch.

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

fprintf(1,['This script will set up the Highly Comparative Time-Series ' ...
                        'Analysis (hctsa) code package from scratch!\n']);
fprintf(1,['We will:' ...
            '\n-1- Add the paths needed for the repository,' ...
            '\n-2- Compile the external time-series toolboxes for this system.\n\n']);

% ------------------------------------------------------------------------------
%% 1. Add the paths:
% ------------------------------------------------------------------------------
fprintf(1,'-1- Adding paths needed for the repository...\n');
input('<<Press any key to continue>>')
try
	startup
catch emsg
	fprintf(1,'error.\n%s\n',emsg.message);
end
fprintf(1,'\n');

% ------------------------------------------------------------------------------
%% Attempt to compile the executables required by the periphery Toolboxes:
% ------------------------------------------------------------------------------
fprintf(1,['\n-2- Compile the binary executables needed for evaluating ' ...
                                                'some operations.\n']);
fprintf(1,['Please make sure that mex is set up with the right compilers for' ...
                                                            ' this system.\n']);
fprintf(1,['Note that errors here are not the end of the world,\nbut mean that ' ...
                        'some operations may fail to execute correctly...\n']);
input('<<Press any key to continue>>')
cd Toolboxes
compile_mex
cd('../');

%-------------------------------------------------------------------------------
%
fprintf(1,'Hope everything compiled ok?!\n\n');

fprintf(1,['All done! Ready when you are to initiate hctsa analysis\nusing a time-series dataset: ' ...
                            'e.g.: TS_init(''INP_test_ts.mat'');\n']);

% Attempt to add a time series
% SQL_add('ts','INP_test_ts.txt')
