function gitInfo = TS_AddGitInfo(whatDir,whatData)
% TS_AddGitInfo   Appends the git information to the hctsa file
%
% Note that after computing, functions and results may change and thus not
% reflect the commit at the time of initailizing (or whenver this is run).
%
%---INPUTS:
% whatDir: the directory to look for a git repository in
% whatData: the data file to append
%
%---OUTPUTS:
% gitInfo, a structure with fields: branch (the current branch),
%           hash (the SHA code for the current commit)
%           remote (the name of the remote)
%           url (the url of the remote)
%
% Also writes output into HCTSA.mat (if whatData is provided)
%
% Uses the getGitInfo function by Andrew Leifer (copyright 2011)

%-------------------------------------------------------------------------------
% Check inputs
%-------------------------------------------------------------------------------
if nargin < 1
    % look in the hctsa directory (one level higher than where TS_compute is)
    calc_dir = which('TS_compute');
    iHere = regexp(calc_dir,'Calculation');
    whatDir = calc_dir(1:iHere-1);
    % whatDir = '.'; % now looky here
end

% Get the gitInfo for the given directory under git version control
gitInfo = getGitInfo();

% Also append this information to the hctsa .mat file provided (if )
if nargin == 2 && ~isempty(whatData)
    % Append to an hctsa data file
    if exist(whatData,'file')
        fileSave = which(whatData);
        save(fileSave,'gitInfo','-append')
        fprintf(1,'Saved git info to %s\n',whatData);
    else
        error('%s does not exist',whatData);
    end
end


%-------------------------------------------------------------------------------
function gitInfo = getGitInfo()
    % Get information about the Git repository in the current directory, including:
    %          - branch name of the current Git Repo
    %          -Git SHA1 HASH of the most recent commit
    %          -url of corresponding remote repository, if one exists
    %
    % The function first checks to see if a .git/ directory is present. If so it
    % reads the .git/HEAD file to identify the branch name and then it looks up
    % the corresponding commit.
    %
    % It then reads the .git/config file to find out the url of the
    % corresponding remote repository. This is all stored in a gitInfo struct.
    %
    % Note this uses only file information, it makes no external program
    % calls at all.
    %
    % This function must be in the base directory of the git repository
    %
    % Released under a BSD open source license. Based on a concept by Marc
    % Gershow.
    %
    % Andrew Leifer
    % Harvard University
    % Program in Biophysics, Center for Brain Science,
    % and Department of Physics
    % leifer@fas.harvard.edu
    % http://www.andrewleifer.com
    % 12 September 2011

    % Copyright 2011 Andrew Leifer. All rights reserved.
    %
    % Redistribution and use in source and binary forms, with or without modification, are
    % permitted provided that the following conditions are met:
    %
    %    1. Redistributions of source code must retain the above copyright notice, this list of
    %       conditions and the following disclaimer.
    %
    %    2. Redistributions in binary form must reproduce the above copyright notice, this list
    %       of conditions and the following disclaimer in the documentation and/or other materials
    %       provided with the distribution.
    %
    % THIS SOFTWARE IS PROVIDED BY <COPYRIGHT HOLDER> ''AS IS'' AND ANY EXPRESS OR IMPLIED
    % WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
    % FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> OR
    % CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    % CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    % SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
    % ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
    % NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
    % ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    %
    % The views and conclusions contained in the software and documentation are those of the
    % authors and should not be interpreted as representing official policies, either expressed
    % or implied, of <copyright holder>.

    gitInfo = struct();
    if ~exist(fullfile(whatDir,'.git'),'file') || ~exist(fullfile(whatDir,'.git','HEAD'),'file')
        % Git is not present
        % warning('No git repository found in %s',whatDir)
        return
    end

    % Read in the HEAD information, this will tell us the location of the file
    % containing the SHA1
    text = fileread(fullfile(whatDir,'.git','HEAD'));
    parsed = textscan(text,'%s');

    if ~strcmp(parsed{1}{1},'ref:') || ~length(parsed{1})>1
        % the HEAD is not in the expected format.
        warning('git HEAD in %s in unexpected format',whatDir)
        return
    end

    path = parsed{1}{2};
    [pathstr, name, ext] = fileparts(path);
    branchName = name;

    % Save branchname
    gitInfo.branch = branchName;

    % Read in SHA1
    SHA1text = fileread(fullfile(whatDir,'.git',pathstr,[name ext]));
    SHA1 = textscan(SHA1text,'%s');
    gitInfo.hash = SHA1{1}{1};

    % Read in config file
    config = fileread(fullfile(whatDir,'.git','config'));
    % Find everything space delimited
    temp = textscan(config,'%s','delimiter','\n');
    lines = temp{1};

    %---------------------------------------------------------------------------
    % Lets find the name of the remote corresponding to our branchName
    remote = '';
    for k = 1:length(lines)

        % Are we at the section describing our branch?
        if strcmp(lines{k},['[branch "' branchName '"]'])
            m = k + 1;
            %While we haven't run out of lines
            %And while we haven't run into another section (which starts with
            % an open bracket)
            while (m<=length(lines) && ~strcmp(lines{m}(1),'[') )
                temp = textscan(lines{m},'%s');
                if length(temp{1})>=3
                    if strcmp(temp{1}{1},'remote') && strcmp(temp{1}{2},'=')
                        %This is the line that tells us the name of the remote
                        remote=temp{1}{3};
                    end
                end
                m = m+1;
            end
        end
    end
    gitInfo.remote = remote;

    %---------------------------------------------------------------------------
    % Find the remote's url
    url = '';
    for k = 1:length(lines)

        %Are we at the section describing our branch?
        if strcmp(lines{k},['[remote "' remote '"]'])
            m = k+1;
            % While we haven't run out of lines
            % And while we haven't run into another section (which starts with
            % an open bracket)
            while (m<=length(lines) && ~strcmp(lines{m}(1),'[') )
                temp=textscan(lines{m},'%s');
                if length(temp{1})>=3
                    if strcmp(temp{1}{1},'url') && strcmp(temp{1}{2},'=')
                        %This is the line that tells us the name of the remote
                        url=temp{1}{3};
                    end
                end
                m = m + 1;
            end
        end
    end
    gitInfo.url = url;

end

end
