% INSTALL_DATABASE   Sets up a mySQL database to store computation results from hctsa

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

fprintf(1,['This script will set up a mySQL database for hctsa (aiding distributed hctsa computation).\n']);
fprintf(1,['We will:\n-1- Set up a connection to a mySQL database,' ...
        '\n-2- Add a default set of time-series analysis operations to the database,']);

% ------------------------------------------------------------------------------
%% Set up the database:
% ------------------------------------------------------------------------------
fprintf(1,'--- Setting up a link to a mySQL database...\n');
reply = input('Do you want to set up a mySQL database to use with hctsa? [''y'' for yes, default: no]','s');
if strcmp(reply,'y')
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
                                ' to a mySQL server to do this\n']);
        % Walks the user through creating the database from a root account and sets
        % up a user account and password:
        SQL_create_db;
        fprintf(1,['Note that the database access ' ...
                    'settings are stored in ''sql_settings.conf''\n']);
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
    fprintf(1,'-2- Populating the database with operations...\n');
    input('<<Press any key to continue>>')

    fprintf(1,'Adding Master operations...\n');
    mop_timer = tic;
    SQL_add('mops','INP_mops.txt','',0)
    fprintf(1,'Master operations added in %s.\n',BF_thetime(toc(mop_timer)));
    clear mop_timer

    fprintf(1,'Adding the default set of operations...\n');
    op_timer = tic;
    SQL_add('ops','INP_ops.txt','',0)

    fprintf(1,'Operations added in %s.\n\n',BF_thetime(toc(op_timer)));
    clear op_timer

    fprintf(1,['All done! Ready when you are to add time series to the database using ' ...
                            'SQL_add(''ts'',''<timeSeries_inputFile.txt>'')...!\n']);
else
    fprintf(1,'\n--- No link between Matlab and a mySQL database will be set up.\n');
    fprintf(1,'You can use TS_init to run hctsa analysis locally within Matlab.\n');
end
