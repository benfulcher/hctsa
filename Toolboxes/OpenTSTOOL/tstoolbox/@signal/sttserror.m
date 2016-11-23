function rs = sttserror(s1, s2)

%tstoolbox/@signal/sttserror
%   Syntax:
%     * rs = sttserror(s1, s2)
%
%   Input Arguments:
%     * s1 - original signal
%     * s2 - predicted signal
%
%   compute error function for prediction of spatial-temporal systems
%   see U. Parlitz "Nonlinear Time Series Analysis", Section 1.10.2.2 Eq.
%   1.10
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,2);

if ~isa(s2,'signal')
        error('s2 must be a signal');
end

if dlens(s1) ~= dlens(s2)
	error('Signals must have same size');
end

s = data(s1);
p = data(s2);

N = dlens(s1, 1);	% Temporal steps
M = dlens(s1, 2);	% Spatial steps

E = zeros(N,1);

for n=1:N
	tau = n:N;
	E = E + [sum((p(tau, :)-s(tau, :)).^2, 2) ./ sum( (s(tau, :)).^2, 2) ; zeros(n-1,1)]; 
end

E  = E ./ (N:-1:1)';

rs = signal(core(E), s1);   % special constructor calling syntax for working routines
a = getaxis(s1, 1);
a = setfirst(a, 0);
rs = setaxis(rs, 1, a);
rs = addhistory(rs, 'Prediction error for spatial systems');
rs = addcommandlines(rs, 's = sttserror(s, s2');
