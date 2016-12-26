function [xpos, u] = firstzero(s)

%tstoolbox/@signal/firstzero
%   Syntax:
%     * [xpos, unit] = firstzero(s)
%
%   Give information about first zero of scalar signal s, using linear
%   interpolation.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt
narginchk(1,1);

if ndim(s) ~= 1
	help(mfilename)
	return
end

d = data(s);
ze = find(d==0);
if isempty(ze)		% no exact zero found, looking for zero crossings
	tmp = d(2:end) .* d(1:end-1);
	ze = find(tmp < 0);
	if isempty(ze)
		error('no zero crossing found')
	else
		A = getaxis(s, 1);
		a = abs(d(ze(1)));
		b = abs(d(ze(1)+1));
		if (a+b) > 0
			xpos = first(A) + (ze(1)-1 + a/(a+b) )*delta(A);
		else
			xpos = first(A) + (ze(1)-1)*delta(A);
		end
		u = unit(A);
	% 	disp(['first zero/zero crossing at x = ' num2str(xpos) ' ' label(A)]);
	end
else		% use first sample that matches zero exactly
	A = getaxis(s, 1);
	xpos = first(A) + (ze(1)-1)*delta(A);
	u = unit(A);
end
