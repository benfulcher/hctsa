function out = benrighttime(tsec)
% Converts the input, tsec, a duration of time in seconds, into an appropriate 
% string for output (i.e., converts to minutes or hours or days as appropriate)
% output is something like '25.5 minutes' or '3.2 days' -- always displays to one decimal
% place

if tsec<=60 % less than a minute, display in seconds
	out = [num2str(round(tsec*10)/10) ' seconds'];
elseif tsec<=60*60 % less than an hour, display in minutes
	out = [num2str(round(tsec/60*10)/10) ' minutes'];
elseif tsec<=60*24*60 % less than a day, display in hours
	out = [num2str(round(tsec/60/60*10)/10) ' hours'];	
% elseif tsec<=60*24*7*60 % less than a week, display in days
else % display in days
	out = [num2str(round(tsec/60/60/24*10)/10) ' days'];
end

end