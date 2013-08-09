% RA_sqlescapestring
% 
% Converts an input string into one that won't interfere with MySQL queries by
% making the required conversions.
% 
% INPUT:
% s, the input string
% 
% OUTPUT:
% S, the output string
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013, Romesh Abeysuriya <romesh.abey@gmail.com>
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function S = RA_sqlescapestring(s)
% Romesh Abeysuriya, 10/12/2012

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

S = regexprep(s,'([\''\\])','\\$1'); % Replace backslashes and single quotes
%S2 = regexprep(s,'([\''\\\_\%])','\\$1'); % Replace backslashes and single quotes

end