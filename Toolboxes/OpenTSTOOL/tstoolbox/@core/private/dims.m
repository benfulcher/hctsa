function d = dims(data)

% Replacement for size 
% Returns only one value for a column vector
% Cuts of any trailing ones

if isempty(data)
	d = [];
else
	d = size(data);
	m = max(find(d ~= 1));

	if ~isempty(m)
		d(max(find(d~=1))+1:end) = [];
	else
		d = [1];
	end
end
