function varargout=first(x)
%FIRST   First three moments.
%    [M S V]=FIRST(X) computes the first three moments mean  M, 
%    standard deviation S and variance V of the multi-column vector X.
%
%    FIRST(...) without output arguments presents the result in
%    a message box.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 1998-2008
% Christian Hoennicke/ Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:32:09 $
% $Revision: 2.2 $
%
% $Log: first.m,v $
% Revision 2.2  2009/03/24 08:32:09  marwan
% copyright address changed
%
% Revision 2.1  2004/11/10 07:07:51  marwan
% initial import
%

error(nargchk(0,3,nargin));
if nargout>3, error('Too many output arguments'), end

for l=1:size(x,2),
    f=find(~isnan(x(:,l)));
    a(l)=mean(x(f,l));
    b(l)=std(x(f,l));
    c(l)=var(x(f,l));
    k(l)=1.96*b(l)/sqrt(length(x));
end     

 
%Ausgabe

if nargout==0
   T=[ ...
   {'Mean:                  ' sprintf('%0.4f ',a)} {''}...
   {'Konfidence Intervalle: ' sprintf('%0.4f  ',k)} {''}...
   {'Standard Deviation:    ' sprintf('%0.4f  ',b)} {''}...
   {'Variance:              ' sprintf('%0.4f ',c)}];
   msgbox(T,'First Moments');
elseif nargout==1
   varargout(1)={a};
elseif nargout==2
   varargout(1)={a};
   varargout(2)={b};
elseif nargout==3
   varargout(1)={a};
   varargout(2)={b};
   varargout(3)={c};
end
