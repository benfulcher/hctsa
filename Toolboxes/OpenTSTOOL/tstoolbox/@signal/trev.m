function rs=trev(s,tau,ntests,method)

%tstoolbox/@signal/trev
%   Syntax:
%     * rs = trev(s,tau,n,method)
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
%       the T[REV] value for the original data set s. The following n
%       values are the T[REV] values for the generated surrogates. There
%       exist a special plothint ('surrbar') for the view function to show
%       this kind of result in the common way.
%
%   This function calculates a special value for the original data set and
%   the n generated surrogate data sets. The T[REV] value is defined as
%   followed:
%
%               [missing equation, see html/pdf documentation]
%
%   In terms of surrogate data test this is a test statistics for time
%   reversibility. The original trev function is located under
%   utils/trev.m and use simple matlab vectors.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

%x1=data(cut(s,1,1+tau,length(data(s))))';


x = [trev(data(cut(s,1,1,length(data(s))))',tau)];

  
for tests=1:ntests
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
  x=[x trev(data(cut(a,1,1,length(data(a))))',tau)];		
end
  

c = core(x');
a=achse(unit(''),0,1);
a=setname(a,'Tests');
rs = signal(c, s);	% special constructor calling syntax for working routines
rs=setyname(rs,'Trev');
rs=setaxis(rs,1,a);
rs = addhistory(rs, ['Computed ' num2str(ntests) ' surrogate data trev values']);
rs = setplothint(rs, 'surrbar');
rs = addcommandlines(rs, 's = trev(s', tau,ntests,method);

return
 