% ------------------------------------------------------------------------------
% install_jconnector
% ------------------------------------------------------------------------------
% 
% Installs the j connector to the current version of Matlab and adds it to the
% java class path
% 
%---INPUT:
% a path specifying a mysql-connector-java file
% 
%---EXAMPLE USAGE:
%
% install_jconnector('Database/mysql-connector-java-5.1.35-bin.jar')
%
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

function install_jconnector(jConnectorWhere)

if nargin < 1
    maybeHere_1 = fullfile('mysql-connector-java-5.1.35-bin.jar');
    maybeHere_2 = fullfile('Database','mysql-connector-java-5.1.35-bin.jar');
    if exist(maybeHere_1)
        jConnectorWhere = maybeHere_1;
    elseif exist(maybeHere_2)
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

jarextDir = fullfile(matlabroot,'java','jarext');

doesExist = exist(fullfile(jarextDir,fileName));

if doesExist
    fprintf(1,'%s already exists in %s.\n\n',fileName,jarextDir);
else
    fprintf(1,'Copying %s to the %s directory...',jConnectorWhere,jarextDir);

    copyfile(jConnectorWhere,jarextDir);

    fprintf(1,' Done.\n\n');
end

% ------------------------------------------------------------------------------
% 2. Add reference to it to the javaclasspath.txt
% ------------------------------------------------------------------------------

doesExist = exist(fullfile(prefdir,'javaclasspath.txt'));

if doesExist
    % Already exists -- 
    fprintf(1,'javaclasspath.txt file already exists in %s\n',prefdir);
    fprintf(1,'You must manually add the following line to it (if it doesn''t already exist):\n');
    fprintf(1,'\n$matlabroot%sjava%sjarext%s%s\n',filesep,filesep,filesep,fileName);
    fprintf(1,'\nRun the following code to edit:\n>> edit %s%sjavaclasspath.txt;\n',prefdir,filesep);
    fprintf(1,'\nOnce this is complete, restart Matlab to start using the mySQL java connector.\n');
else
    % Easy, we can create it:
    fprintf(1,'Writing new javaclasspath.txt file in the preference directory: %s ...',prefdir);
    fid = fopen(fullfile(prefdir,'javaclasspath.txt'),'w');
    fprintf(fid,'$matlabroot%sjava%sjarext%s%s\n',filesep,filesep,filesep,fileName);
    % fprintf(fid,'%s%s%s\n',jarextDir,filesep,fileName);
    fclose(fid);
    fprintf(1,' Done.\n Success! Restart Matlab to start using the mySQL java connector.\n');
end

end