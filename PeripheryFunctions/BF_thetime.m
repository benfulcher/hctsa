% BF_thetime
% 
% Converts the input, tsec, a duration of time in seconds, into an appropriate
% string for output (i.e., converts to minutes or hours or days as appropriate)
% output is something like '25.5 minutes' or '3.2 days' -- always displays to
% one decimal place.
% 
% INPUTS:
% tsec, the duration in seconds
% formatlong, (i) 0: display short units (like 's' instead of 'seconds')
%                    [default]
%             (ii) 1: display long units (like 'seconds' instead of 's')
% 
% OUTPUT:
% timestring, an interpretable text version of the input time.
% 
% This code is useful for displaying user feedback on tic/toc statements.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
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

function timestring = BF_thetime(tsec,formatlong)
% Ben Fulcher, 2009

if nargin < 2 || isempty(formatlong)
    formatlong = 0; % set to 1 to use a longer format for the unit
end

if tsec < 1E-3
    if formatlong
        timestring = '< 1 milliseconds';
    else
        timestring = '< 1ms';
    end
elseif tsec < 1 % less than a second, display in integer number of milliseconds
    if formatlong
        timestring = sprintf('%.0f milliseconds',tsec*1000);
    else
        timestring = sprintf('%.0fms',tsec*1000);
    end
elseif tsec <= 60 % less than a minute, display in seconds
    if formatlong
        timestring = sprintf('%.1f seconds',tsec);
    else
    	timestring = sprintf('%.1fs',tsec);
    end
elseif tsec <= 60*60 % less than an hour, display in minutes
    if formatlong
        timestring = sprintf('%.1f minutes',tsec/60);
    else
        timestring = sprintf('%.1fmin',tsec/60);
    end
elseif tsec <= 60*24*60 % less than a day, display in hours
    if formatlong
        timestring = sprintf('%.1f hours',tsec/60/60);
    else
        timestring = sprintf('%.1fh',tsec/60/60);
    end
% elseif tsec<=60*24*7*60 % less than a week, display in days
else % display in days
	timestring = sprintf('%.1f days',tsec/60/60/24);
end

end