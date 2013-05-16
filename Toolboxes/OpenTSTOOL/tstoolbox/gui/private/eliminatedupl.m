function ret = eliminatedupl(data)

% Elimiert doppelt vorkommende Strings 
% einer Liste (als Cell-Array)

ret = {};

if ~isempty(data) 
	n = length(data);	
	work =  data{1};
	ret{1} = work;
	for i=2:n
		if ~strcmp(data{i}, work)	
			work = 	data{i};
			ret{end+1} = work;
		end
	end
end






