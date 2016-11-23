function rs = plosivity(s, blen, flen, thresh, windowtype)

%tstoolbox/@signal/plosivity
%   Syntax:
%     * rs = plosivity(s, blen) => flen=1 , thresh=0, windowtype = 'Rect'
%     * rs = plosivity(s, blen, flen) => thresh=0, windowtype = 'Rect'
%     * rs = plosivity(s, blen, flen, thresh) => windowtype = 'Rect'
%     * rs = plosivity(s, blen, flen, thresh, windowtype)
%
%   Compute plosivity of a spectrogram. See also: window for list of
%   possible window types.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,5)

if nargin < 3
	flen = 1;
end
if nargin < 4
	thresh = 0;
end
if nargin < 5
	windowtype =  'Rect';
end

bweight = window(blen, windowtype);
fweight = window(flen, windowtype);

rs = signal(core(plosivity(data(s),bweight, fweight, thresh)), s);	
rs = setaxis(rs, 1, getaxis(s, 2));
rs = addhistory(rs,  ['Computed Plosivity']);
rs = addcommandlines(rs, 's = plosivity(s', blen, flen, thresh, windowtype);

