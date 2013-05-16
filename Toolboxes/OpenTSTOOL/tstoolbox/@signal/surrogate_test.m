function rs=surrogate_test(s1,ntests,method,func)

%tstoolbox/@signal/surrogate_test
%   Syntax:
%     * rs=surrogate_test(s, ntests, method,func)
%
%   Input Arguments:
%     * s - has to be a real, scalar signal
%     * ntests - is the number of surrogate data sets to create
%     * method - method to generate surrogate data sets:
%          + 1: surrogate1
%          + 2: surrogate2
%          + 3: surrogate3
%     * func - string with matlab-code, have to return a signal object
%       with a scalar time series. The data to process is a signal object
%       referred by the qualifier s (see example).
%
%   Output Arguments:
%     * rs is a signal object with a three dimensional time series. The
%       first component is the result of the func function applied to the
%       original data set s. The second component is the mean of the
%       result of the func function applied to the ntests surrogate data
%       sets. The third component is the standard deviation. There is a
%       special plothint ('surrerrorbar') for the view function to show
%       this result in the common way.
%
%   surrogate_test runs an automatic surrogate data test task. It
%   generates ntests surrogate data sets an performs the func function to
%   each set. func is a string with matlab-code who returns a signal s
%   with a scalar time series.
%
%   Example:
%st = surrogate_test(s, 10, 1, 1, 'largelyap(embed(s,3,1,1), 128,20,10);');
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt



s=s1;
s=eval(func);
x=[data(s)];
for tests=1:ntests
  switch(method)
   case 1
    s=surrogate1(s1);
   case 2
    s=surrogate2(s1);
   case 3
    s=surrogate3(s1);
   otherwise
    s=surrogate1(s1);
  end
  s=eval(func);
  x=[x data(s)];
end

a=[];
for test=1:length(x(:,1))
  a=[a; mean(x(test,2:end)) std(x(test,2:end))];
end

a=[x(:,1) a];



rs = signal(core(a),s1);	% special constructor calling syntax for working routines
%rs=setaxis(rs,1,a);
rs = addhistory(rs, ['Computed ' num2str(ntests) ' surrogate data function values']);
rs = setplothint(rs, 'surrerrorbar');
rs = addcommandlines(rs, 's = surrogate_test(s', ntests,method,func);

return

