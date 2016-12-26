function install_jconnector(jConnectorWhere,permanentDir)
% install_jconnector     Installs the j connector to the current version of
% Matlab and adds it to the java class path
%
%---INPUTS:
% jConnectorWhere: a path specifying a mysql-connector-java file
% permanentDir: the directory to install the mysql java connector to (by default
%               this will be the jarext directory of the current Matlab installation)
%               Consider specifying this to be a custom directory of java files
%               (with a static path) if you don't have permission to install in
%               the jarext directory.
%
%---EXAMPLE USAGE:
%
% Installs the j-connector in the jarext directory of the Matlab installation:
% install_jconnector('Database/mysql-connector-java-5.1.35-bin.jar')

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

if nargin < 1 || isempty(jConnectorWhere)
    maybeHere_1 = fullfile(pwd,'mysql-connector-java-5.1.35-bin.jar');
    maybeHere_2 = fullfile(pwd,'Database','mysql-connector-java-5.1.35-bin.jar');
    if exist(maybeHere_1,'file')
        jConnectorWhere = maybeHere_1;
    elseif exist(maybeHere_2,'file')
        jConnectorWhere = maybeHere_2;
    else
        error('Please specify the location of the mysql-connector-java file');
    end
end

fprintf(1,'Installing mySQL j-connector to allow Matlab to communicate with a mySQL server.\n');

fileName = regexp(jConnectorWhere,filesep,'split');
fileName = fileName{end};

% ------------------------------------------------------------------------------
% 1. Copy the mysql-connector to the jar/jarext directory in matlabroot
% ------------------------------------------------------------------------------

if nargin < 2 || isempty(permanentDir)
    permanentDir = fullfile(matlabroot,'java','jarext'); % the jarext directory
    fprintf(1,'Installing the j-connector to %s by default\n',permanentDir);
else
    % Check the provided directory exists
    if ~exist(permanentDir,'dir')
        error('%s is not a valid target directory',permanentDir);
    end
end

doesExist = exist(fullfile(permanentDir,fileName),'file');

if doesExist
    fprintf(1,'%s already exists in %s.\n\n',fileName,permanentDir);
else
    input(sprintf('Copying %s\nto the %s directory...\n[press any key to continue]',jConnectorWhere,permanentDir),'s');

    try
        copyfile(jConnectorWhere,permanentDir);
    catch emsg
        fprintf(1,'\n');
        error(['Error copying file %s to %s.\n(Write permission error?).\nYou should either provide a valid custom ' ...
                    'directory (to which you have write access),' ...
                    '\nor do this step manually (see Documentation).'],...
                    jConnectorWhere,permanentDir);
    end

    fprintf(1,' Done.\n\n');
end

% ------------------------------------------------------------------------------
% 2. Add reference to it to the javaclasspath.txt
% ------------------------------------------------------------------------------

doesExist = exist(fullfile(prefdir,'javaclasspath.txt'),'file');

if doesExist
    % Already exists --
    fprintf(1,'javaclasspath.txt file already exists in %s\n',prefdir);
    fprintf(1,'You must manually add the following line to it (if it doesn''t already exist):\n');
    if nargin < 2
        fprintf(1,'\n$matlabroot%sjava%sjarext%s%s\n',filesep,filesep,filesep,fileName);
    else
        fprintf(1,'\n%s%s%s\n',permanentDir,filesep,fileName);
    end
    fprintf(1,'\nRun the following code to edit:\n>> edit %s%sjavaclasspath.txt;\n',prefdir,filesep);
    fprintf(1,'\nOnce this is complete, restart Matlab to start using the mySQL java connector.\n');
else
    % Easy, we can create it:
    fprintf(1,'Writing new javaclasspath.txt file in the preference directory: %s ...',prefdir);
    fid = fopen(fullfile(prefdir,'javaclasspath.txt'),'w');
    if nargin < 2
        fprintf(fid,'$matlabroot%sjava%sjarext%s%s\n',filesep,filesep,filesep,fileName);
    else
        fprintf(fid,'\n%s%s%s\n',permanentDir,filesep,fileName);
    end
    % fprintf(fid,'%s%s%s\n',permanentDir,filesep,fileName);
    fclose(fid);
    fprintf(1,' Done.\nSuccess! Restart Matlab to start using the mySQL java connector.\n');
end

end
