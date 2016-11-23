function s = addaxis(s, a)

% signal/addaxis
% Add one new axis to signal object
%
% s = addaxis(s)    		% add default axis
% s = addaxis(s, axis)		% add axis given by argument

narginchk(1,2);

if nargin < 2
	s.xaxes{end+1} = achse;
else
	s.xaxes{end+1} = a;
end
