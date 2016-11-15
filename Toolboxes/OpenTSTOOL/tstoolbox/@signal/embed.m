function rs = embed(s, dim, delay, shift, windowtype)

%tstoolbox/@signal/embed
%   Syntax:
%     * emb = embed(s, dim, delay, shift, windowtype)
%
%   Input arguments:
%     * dim - embedding dimension
%     * delay - time delay (optional)
%     * shift - shift for two sequent time delay vectors (optional)
%     * windowtype - type of window (optional)
%
%   Output arguments:
%     * emb - n by dim array, each row contains the coordinates of one
%       point
%
%   Embeds signal s with embedding dimension dim and delay delay (in
%   samples). s must be a scalar time series. The default values for dim
%   and delay are equal to one. The default value for windowtype is
%   'Rect', which is currently the only possible value.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2, 5);

if nargin < 3, delay = 1; end
if nargin < 4, shift = 1; end
if nargin < 5, windowtype = 'Rect'; end

if dlens(s,2) > 1
	error('Cannot embed multi-dimensional time series');
end

c = embed(s.core, dim, delay,shift, windowtype); 	% call real working routine for parent core object
rs = signal(c, s);				% special constructor calling syntax for working routines

a = getaxis(s,1);
rs = setaxis(rs,1,setdelta(a,delta(a)*shift));
a = achse(unit,1,1);
a = setname(a,'Embedding dimension');
rs = setaxis(rs,2,a);
rs = addhistory(rs,'Time delay reconstruction');
rs = addcommandlines(rs,'s = embed(s',dim,delay,shift,windowtype);

switch dim
	case 2
		rs = setplothint(rs, 'xyplot');
	case 3
		rs = setplothint(rs, '3dcurve');
end

end
