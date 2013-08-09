% RA_keyboard
% 
% Romesh Abeysuriya's replacement of Matlab's 'keyboard' command
% Keyboard debug caller 
% Provides more information including the stack trace than simply using 
% 'keyboard'
% Use RA_keyboard() directly instead of keyboard()
% Note this function will not work if it is the last line in the program
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

function RA_keyboard
% By Romesh Abeysuriya, 15-11-12

[stack_trace] = dbstack;
fprintf(1,'\n')
for j = 2:length(stack_trace)
    fprintf('%s:%i << ',stack_trace(j).name,stack_trace(j).line);
end
fprintf('\b\b\b\b\n');
eval(sprintf('dbstop in %s at %u',stack_trace(2).name,stack_trace(2).line+1))

end