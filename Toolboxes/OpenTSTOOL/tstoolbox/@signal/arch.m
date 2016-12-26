function [rs, archetypes] = arch(s, na, mode)

%tstoolbox/@signal/arch
%   Syntax:
%     * [rs, archetypes]=arch(s, na, mode='normalized')
%
%   Input arguments:
%     * na - number of generated archetypes
%     * mode - mode can be one of the following : 'normalized' , 'mean',
%       'raw' (optional)
%
%   Archetypal analysis of column orientated data set:
%     * each row of data is one 'observation', e.g. the sample values of
%       all channels in a multichannel measurement at one point in time
%     * in mode 'normalized' each column of data is centered by removing
%       its mean and then normalized by dividing through its standard
%       deviation before the covariance matrix is calculated
%     * in mode 'mean' only the mean of every column of data is removed
%     * in mode 'raw' no preprocessing is applied to data
%
%   Default value for mode is 'normalized'.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,3);

if nargin < 3
	mode = 'normalized';
end

if ndim(s)~=2
	error('arch needs a signal with two dimensions as input')
end

dat = data(s);
n = dlens(s,1);
m = dlens(s,2);

htext = {''};
htext{end+1} = 'Applied Archetypal Analysis (arch)';

[Z,A,B] = arch2(data(s)', na, mode)

a = achse(unit, 1,1);
a = setname(a, 'Archetyp');
 
rs = signal(core(A'), s);
rs = addhistory(rs, htext);
rs = addcommandlines(rs, 's = arch(s', na, mode);
rs = setaxis(rs, 2, a);
rs = setplothint(rs, 'multigraph');

archetypes = signal(core(Z),s);
archetypes = addhistory(archetypes, htext);
archetypes = addcommandlines(archetypes, '[dummy, s] = arch(s', na, mode);
archetypes = setaxis(archetypes, 2, a);
archetypes = setaxis(archetypes, 1, getaxis(s,2));
