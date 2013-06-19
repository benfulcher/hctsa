function timestring = benrighttime(tsec)
% Converts the input, tsec, a duration of time in seconds, into an appropriate 
% string for output (i.e., converts to minutes or hours or days as appropriate)
% output is something like '25.5 minutes' or '3.2 days' -- always displays to one decimal
% place
% Ben Fulcher, 2009

if tsec < 1E-3
    timestring = '< 1 millisecond';
elseif tsec < 1 % less than a second, display in integer number of milliseconds
    timestring = sprintf('%.0f milliseconds',tsec*1000);
elseif tsec <= 60 % less than a minute, display in seconds
	timestring = sprintf('%.1f seconds',tsec);
elseif tsec <= 60*60 % less than an hour, display in minutes
    timestring = sprintf('%.1f minutes',tsec/60);
elseif tsec <= 60*24*60 % less than a day, display in hours
	timestring = sprintf('%.1f hours',tsec/60/60);
% elseif tsec<=60*24*7*60 % less than a week, display in days
else % display in days
	timestring = [num2str(round(tsec/60/60/24*10)/10) ' days'];
end

end