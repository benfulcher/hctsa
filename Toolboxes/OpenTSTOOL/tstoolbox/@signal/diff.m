function rs = diff(s, nth)

%tstoolbox/@signal/diff
%   Syntax:
%     * diff(s, nth)
%
%   Compute the nth numerical derivative along dimension 1. s has be to
%   sampled equidistantly.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,2);

if nargin < 2
	nth = 1;
end

a = getaxis(s,1);
switch  resolution(a)
	case 'linear'
		c = diff(s.core, nth, delta(a)); 		% call real working routine for parent core object
		rs = signal(c, s);				% special constructor calling syntax for working routines
		rs = setyunit(rs, yunit(s)/(unit(a)^nth));
		a = setfirst(a, first(a) + nth*delta(a)/2);
		rs = setaxis(rs,  1, a);
		rs = addhistory(rs,  ['Computed ' num2str(nth) '# numerical derivative along dimension 1'] );
		rs = addcommandlines(rs, 's = diff(s', nth);
	
        otherwise
	   x=spacing(s)';
	   x2=x(1:end-1);
	   x=x(2:end);
	   
	   y=data(s);
	   y2=y(1:end-1,:);
	   y=y(2:end,:);
	   y=(y-y2);
	   for i=1:length(y(1,:))
	     y(:,i)=y(:,i)./(x-x2);
	   end
	   x=(x-x2)/2+x2;
	   
	   rs=signal(core(y),s);
	   rs=setplothint(rs,'multigraph');
	   rs = setyname(rs, 'd ld N(r)/d lg r');
	   rs = setyunit(rs, yunit(s)/(unit(a)^nth));
	   rs = setaxis(rs,  1, achse(x));
	   rs = addhistory(rs,  ['Computed ' num2str(nth) '# numerical derivative along dimension 1'] );
	   rs = addcommandlines(rs, 's = diff(s', nth);
	   
	   
%	   error('Data values are not sampled equidistantly');
	 
end


