function d = setlabel(d, n)

%tstoolbox/@description/setlabel
%   Syntax:
%     * d = setlabel(d, label)
%
%   the label field of a description is used to give a signal some 'tag'
%   which remains constant through various processing steps
%   e.g. which topic this signal belongs to
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,2);

d.label = n;

