function rs = delaytime(s, maxdelay, past)

%tstoolbox/@signal/delaytime
%   Syntax:
%     * tau = delaytime(s, maxdelay, past)
%
%   Input arguments:
%     * maxdelay - maximal delay time
%     * past - ?
%
%   Compute optimal delaytime for a scalar timeseries with method of
%   Parlitz and Wichard.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(3,3);

if ndim(s) ~= 1
	help(mfilename)
	return
end


ITERATIONS = 64;


len = dlens(s, 1)-maxdelay;
x = data(s);
[y, index] = sort(x(1:len));

err = zeros(maxdelay+1, 1);

for i=1:ITERATIONS
	while(1)
		ref = ceil(rand(1,1)*len);
		actual = index(ref);
    	pre = index(find(abs(index(1:ref-1)-actual) > past));
		post = index(find(abs(index(ref+1:end)-actual) > past));
		if ~isempty(pre) & ~isempty(post)
			pre = pre(end);
			post = post(end);		
			break
		end
	end
	%x([pre actual post])

	%errpre = abs((x(pre:pre+maxdelay)-x(actual:actual+maxdelay))/(x(pre)-x(actual)));
	%errpost = abs((x(post:post+maxdelay)-x(actual:actual+maxdelay))/(x(post)-x(actual)));
	errpre = abs(x(pre:pre+maxdelay)-x(actual:actual+maxdelay));
	errpost = abs(x(post:post+maxdelay)-x(actual:actual+maxdelay));
	err = err + errpre  + errpost;
end
 
rs = signal(core(err/ITERATIONS), s);	% special constructor calling syntax for working routines
a = getaxis(s, 1); 
dl = delta(a);
a = setfirst(a, 0);
rs = setaxis(rs, 1, a);
rs = setyunit(rs, unit);		% acf values are scalars without unit
rs = addhistory(rs, 'Delaytime plot');
rs = addcommandlines(rs, 's = delaytime(s', maxdelay, past);
