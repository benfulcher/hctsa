function MakeTableOfFunctions()
% MakeTableOfFunctions Outputs a csv summarizing INP_mops.txt
% ------------------------------------------------------------------------------
% Copyright (C) 2017, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems (2017).
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

% Read in INP_mops.txt to get a list of functions
fid = fopen('INP_mops.txt');
Mops = textscan(fid,'%s %s','CommentStyle','#','CollectOutput',1);
fclose(fid);
Mops = Mops{1};
Mops = Mops(:,1);

% Split at the parenthesis:
MopsSplit = regexp(Mops,'(','split','ignorecase');
mFiles = cellfun(@(x)x{1},MopsSplit,'UniformOutput',false);
mFiles = unique(mFiles);
numMFiles = length(mFiles);

% Get help lines as descriptions:
helpLines = cell(numMFiles,1);
for i = 1:numMFiles
    helpMe = help(mFiles{i});
    newLines = regexp(helpMe,'\n');
    theTopLine = helpMe(1:newLines(1)-1);
    notSpaces = regexp(theTopLine,'\w*');
    helpLines{i} = theTopLine(notSpaces(2):end);
end

% Make a table
T = table(mFiles,helpLines);

% Output to csv:
writetable(T,'codeSummary.csv','Delimiter',',','QuoteStrings',true);

end
