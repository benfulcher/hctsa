function s = odesolver(funcname, initcond, length, samplerate, solver, yunit)

% signal/odesolver
% obtain a signal by integrating a set of ordinary differntial equations
%
% si = odesolver(funcname, initcond, length, samplerate, solver, yunit)
%
% funcname is the name of the function to integrate
% initcond is a vector of start values for the integration
% length is measured in seconds
% samplerate in Hertz
% yunit is an optional argument

if nargin < 2
	help(mfilename)
	return	
end
if nargin < 3
	length = 1; % one second
end
if nargin < 4
	samplerate = 8000;	% Hertz
end
if nargin < 5
	solver = 'ode45';
end
if nargin < 6
	yunit = unit;	
end

[T,Y] = feval(solver, funcname, 0:1/samplerate:length, initcond);

s = signal(Y, ['ODE integration with initital conditions ' num2str(initcond(:)')]); % '
s = setaxis(s, 1,achse(unit('s'), 0, 1/samplerate));
s = setyunit(s, yunit);

if dlens(s,2)==2
	s = setplothint(s, 'xyplot');
elseif dlens(s,2)==3
	s = setplothint(s, '3dcurve');
end
