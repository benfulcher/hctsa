function rs=tc3(s,tau,n,method)

%tstoolbox/@signal/tc3
%   Syntax:
%     * rs = tc3(s,tau,n,method)
%
%   Input Arguments:
%     * tau - see explaination below
%     * n - number of surrogate data sets to generate
%     * method - method to generate the surrogate data sets:
%          + 1: surrogate1
%          + 2: surrogate2
%          + 3: surrogate3
%
%   Output Arguments:
%     * rs is a row vector, returned as signal object. The first item is
%       the T[C3] value for the original data set s. The following n
%       values are the T[C3] values for the generated surrogates. There
%       exist a special plothint ('surrbar') for the view function to show
%       this kind of result in the common way.
%
%   This function calculates a special value for the original data set and
%   the n generated surrogate data sets. The T[C3] value is defined as
%   followed:
%
%               [missing equation, see html/pdf documentation]
%
%   In terms of surrogate data test this is a test statistics for higher
%   order moments. The original tc3 function is located under utils/tc3.m
%   and use simple matlab vectors.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

%x1=data(cut(s,1,1+tau,length(data(s))))';

x = [tc3(data(cut(s,1,1+tau,length(data(s))))',tau)];

for i=1:n
  switch(method)
   case 1
    a=surrogate1(s);
   case 2
    a=surrogate2(s);
   case 3
    a=surrogate3(s);
   otherwise
    a=surrogate1(s);
  end
  x=[x tc3(data(cut(a,1,1+tau,length(data(a))))',tau)];		
end

c = core(x');
a=achse(unit(''),0,1);
a=setname(a,'Tests');
rs = signal(c, s);	% special constructor calling syntax for working routines
rs=setyname(rs,'Tc3');
rs=setaxis(rs,1,a);
rs = addhistory(rs, ['Computed ' num2str(n) ' surrogate data tc3 values']);
rs = setplothint(rs, 'surrbar');
rs = addcommandlines(rs, 's = tc3(s', tau,n,method);


return
 