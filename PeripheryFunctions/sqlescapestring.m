function S = sqlescapestring(s)
% Input s, turns into a string that won't interfere with MySQL queries
% Romesh Abeysuriya 10/12/2012

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

S = regexprep(s,'([\''\\])','\\$1'); % Replace backslashes and single quotes
%S2 = regexprep(s,'([\''\\\_\%])','\\$1'); % Replace backslashes and single quotes

