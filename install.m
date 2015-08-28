% INSTALL   Installs the hctsa code package from scratch.

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

fprintf(1,['This script will set up the Highly Comparative Time-Series ' ...
                                'Analysis code package from scratch!\n'])
fprintf(1,['In the following order, we will:' ...
            '\n-1- Add the paths needed for the repository,' ...
            '\n-2- Set up a connection to a mySQL database,' ...
            '\n-3- Add the time-series analysis operations to the database,' ...
            '\n-4- Compile the external time-series toolboxes for this system.\n\n'])

% ------------------------------------------------------------------------------
%% 1. Add the paths:
% ------------------------------------------------------------------------------
fprintf(1,'-1- Adding paths needed for the repository...\n');
input('<<Press any key to continue>>')
try
	startup
catch emsg
	fprintf(1,'error.\n%s\n',emsg.message)
end
fprintf(1,'\n');

% ------------------------------------------------------------------------------
%% 2. Set up the database:
% ------------------------------------------------------------------------------
fprintf(1,'-2- Setting up a link to a mySQL database...\n');
reply = input('Do you want to set up a mySQL database to use with hctsa? [''y'' for yes, default: no]','s');
if strcmp(reply,'y')
    choseDatabase = 1;
    reply = '';
    while isempty(reply)
        reply = input('Do you need help setting up a mySQL database? [y/n]','s');
        if ~ismember(reply,{'y','n'}) && ~isempty(reply)
            fprintf(1,'Unknown response: ''%s''.\n',reply);
            reply = '';
        end
    end
    if strcmp(reply,'y') % Set up mySQL database
        fprintf(1,['Setting up the database now--NB: you need to have root access' ...
                                ' to a mySQL server to do this\n'])
        % Walks the user through creating the database from a root account and sets
        % up a user account and password:
        SQL_create_db;
        fprintf(1,['Note that the database access ' ...
                    'settings are stored in ''sql_settings.conf''\n'])
    else
        fprintf(1,['Ok, then I''ll assume you''ve configured your ''sql_settings.conf'' file correctly.\n' ...
                    '(see Documentation for details on how to configure your Matlab/mySQL connection).\n']);
    end

    % ------------------------------------------------------------------------------
    %% 3. Create all (empty) tables in the database
    % ------------------------------------------------------------------------------
    SQL_create_all_tables;
    fprintf(1,'\n');

    % ------------------------------------------------------------------------------
    %% 4. Populate the new tables with operations
    % ------------------------------------------------------------------------------
    fprintf(1,'-3- Populating the database with operations...\n')
    input('<<Press any key to continue>>')

    fprintf(1,'Adding Master operations...\n');
    mop_timer = tic;
    SQL_add('mops','Database/INP_mops.txt','',0)
    fprintf(1,'Master operations added in %s.\n',BF_thetime(toc(mop_timer)))
    clear mop_timer

    fprintf(1,'Adding all operations...\n');
    op_timer = tic;
    SQL_add('ops','Database/INP_ops.txt','',0)
    % fprintf(1,'Adding a reduced set of operations...\n'); op_timer = tic;
    % SQL_add('ops','Database/INP_ops_reduced.txt','',0)
    fprintf(1,'Operations added in %s.\n\n',BF_thetime(toc(op_timer)))
    clear op_timer
else
    choseDatabase = 0;
    fprintf(1,['\n-3- No link between Matlab and a mySQL database will be set up.\n'])
    fprintf(1,['You will need to use TS_init to run hctsa analysis within Matlab locally.\n'])
    input('<<Press any key to continue>>')
end

% ------------------------------------------------------------------------------
%% 5. Attempt to compile the executables required by the periphery Toolboxes:
% ------------------------------------------------------------------------------
fprintf(1,['\n-4- Attempting to compile the binary executables needed for evaluating ' ...
                                                        'some operations.\n'])
fprintf(1,['Please make sure that mex is set up with the right compilers for' ...
                                                            ' this system.\n'])
fprintf(1,['Note that errors here are not the end of the world,\nbut mean that ' ...
                        'some operations may fail to execute correctly...\n'])
input('<<Press any key to continue>>')
cd Toolboxes
compile_mex
cd('../');

fprintf(1,'Hope everything compiled ok?!\n\n')
if choseDatabase
    fprintf(1,['All done! Ready when you are to add time series to the database using ' ...
                            'SQL_add(''ts'',''<timeSeries_inputFile.txt>'')...!\n']);
else
    fprintf(1,['All done! Ready when you are to initiate hctsa analysis\nusing a time-series dataset: ' ...
                            'e.g.: TS_init(''INP_test_ts.mat'');\n']);
end

% Attempt to add a time series
% SQL_add('ts','INP_test_ts.txt')
