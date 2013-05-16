function a = setunit(a,u)

%tstoolbox/@achse/setunit
%   Syntax:
%     * a = setunit(a,u)
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

a.unit = u;
a.quantity = quantity(u);
a.name = a.quantity;
