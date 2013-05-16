function S = sqlescapestring(s)
% Input s, turns into a string that won't interfere with MySQL queries
% Ben Fulcher 30/11/2009

% From MySQL manual pages:
% http://dev.mysql.com/doc/refman/4.1/en/string-syntax.html
% \'	A single quote (“'”) character -- also ''
% \"	A double quote (“"”) character.
% \b	A backspace character.
% \n	A newline (linefeed) character.
% \r	A carriage return character.
% \t	A tab character.
% \Z	ASCII 26 (Control-Z). See note following the table.
% \\	A backslash (“\”) character.
% \%	A “%” character. See note following the table.
% \_	A “_” character. See note following the table.

% s

% Backslashes (must be done first -- or else picking up on backslashed added by commands)
a = strfind(s,'\');
if ~isempty(a), s = SUB_addbackslash(s,a); end

% Single quotes
a = strfind(s,'''');
if ~isempty(a), s = SUB_addbackslash(s,a); end

% Percentage signs
% a = strfind(s,'%');
% if ~isempty(a), s = SUB_addbackslash(s,a); end

% Underscores
% a = strfind(s,'_');
% if ~isempty(a), s = SUB_addbackslash(s,a); end

S = s;

function s1 = SUB_addbackslash(s1,a)
	a = sort(a,'descend'); % so insertions don't change the string indicies
	for i=1:length(a)
		s1 = [s1(1:a(i)-1) '\' s1(a(i):end)];
	end
end


end