function [c_out,s_out]=xcf(varargin)
%XCF   Computes and plots crosscorrelation.
%    [C,S]=XCF(X,Y [,T,FLAG]) computes crosscorrelation C
%    between the data in the vectors X and Y with the maximal 
%    lag T (optional). If FLAG is set, the plot is suppressed. 
%    The optional output S stores the 5% significance level.
%
%    See also COR, COV, XCORR, ACF

% Copyright (c) 2001-2002 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2007/12/20 16:26:07 $
% $Revision: 2.2 $
%
% $Log: xcf.m,v $
% Revision 2.2  2007/12/20 16:26:07  marwan
% changed gpl splash behaviour
%
% Revision 2.1  2004/11/10 07:09:40  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

norm = 'nonorm';
error(nargchk(1,5,nargin));

i_num=find(cellfun('isclass',varargin,'double'));
x = varargin{i_num(1)};
t = []; flag = 0;
if length(i_num) > 1; y = varargin{i_num(2)}; end
if length(i_num) > 2; t = varargin{i_num(3)}; end
if length(i_num) > 3; flag = varargin{i_num(4)}; end

% check the input
if length(i_num)<2
  y=x;
end
if length(x)~=length(y)
  error('The lengths of x and y must match.')
end


i_char=find(cellfun('isclass',varargin,'char'));
if ~isempty(i_char)
    norm = varargin{i_char};
end
 


% check if xcorr exist
if isempty(which('xcorr'))
  disp('Sorry, signal processing toolbox needed.')
  if nargout==1
    c_out(1:t)=NaN;
  elseif nargout==2
    c_out(1:t)=NaN;
    s_out(1:t)=NaN;  
  end

else

% compute xcf
if strcmpi(norm,'norm')
    x=trafo(x);
    y=trafo(y);
end
x = (x - mean(x)) ./ std(x);
y = (y - mean(y)) ./ std(y);
[c, s]=xcorr(x,y,t,'coeff');
c=flipud(c);
s=s';
L=length(x)-abs(s);
si=2./sqrt(L+2);

% output
if flag==0
bar(s,c,'k')
hold on
plot(s,si,'r-.')
plot(s,-si,'r-.')
line([0 0],[-1 1],'linestyle',':')
ha=axis;
line([ha(1) ha(2)],[0 0],'linestyle',':')
hold off
end

if nargout==1
  c_out=c;
elseif nargout==2
  c_out=c;
  s_out=si;  
end

end
