function s = tssize(data)

% function s = tssize(data)
% Replacement for size 
% Returns only one value for a columns vector
% Cuts of any trailing ones

s = size(data);

while s(end) == 1
	s = s(1:end-1);
end

