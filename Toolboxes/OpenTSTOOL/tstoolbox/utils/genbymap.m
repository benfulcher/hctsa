function s = genbyode(system, len, initcond, params, yunit)

% signal/genbymap
%
% Generate signal by an iterated map,
% 
% si = genbyode(map, len, initcond, params, yunit)
%
% len - length of signal, in seconds
% initcond - vector of inital conditions
% params - vector of parameter values
%
% yunit is an optional argument
%
% system may be one of the following :
%	'Henon'
%	'Baker'
%	'Tentmap'
%
%
% cmerk 1999

if nargin < 1
	system = 'Henon';		
end
if nargin < 2
	len = 500;			
end
if nargin < 3
	initcond = [0.13 0.19];
end
if nargin < 4
	params = [-1.4 0.3];
end
if nargin < 5
	yunit = unit;		% no dimension
end

optinaltext = '';

switch system
	case	'Henon'
		tmp = henon(len, [params initcond]);
		optinaltext = [''];
	case	'Baker'
		tmp = baker(len, [params initcond]);
		optinaltext = [''];
	case	'Tentmap'
		tmp = tentmap(len, [params initcond]);
		optinaltext = [''];				
	otherwise	
		error(['system ' map ' not supported']);
end	
	
s = signal(tmp, ['Generated ' system ' signal' optinaltext]);
s = setyunit(s, yunit);

if (size(tmp, 2) == 2) 
	s = setplothint(s, 'xypoints');
end
if (size(tmp, 2) == 3) 
	s = setplothint(s, '3dpoints');
end

