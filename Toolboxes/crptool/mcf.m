function [a_out, b_out]=mcf(x,y,w,T)
%MCF   Plots maximal correlation function.
%    MCF(X,Y [,W,T]) plots the maximal correlation function up to
%    the maximal lag of T and by using a boxcar window size of 2*W+1.  
%    If W=[], the default boxcar window size is 11.
%
%    MCOR=MCF(...) suppresses the plot and puts the maximal 
%    correlation function into the vector MCOR.
%
%    [TIME, MCOR]=MCF(...) suppresses the plot and puts the maximal 
%    correlation function into the vector MCOR and its time axis 
%    in the vector TIME.
%
%    Examples: x = sin(0:.05:10) + .5*randn(1,201);
%              y = cos(0:.05:10);
%              mcf(x,y,[],20)
%
%    References:
%    Breiman, L., Friedman, J. H.:
%    Estimating Optimal Transformations for Multiple regression
%    and Correlation, J. Am. Stat. Assoc., Vol. 80, No. 391, 
%    1985.
%    Voss, H., Kurths, J.: 
%    Reconstruction of nonlinear time delay models from data by the 
%    use of optimal transformations, Phys. Lett. A, 234, 1997.
%
%    See also: ACE

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 2001-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:32:57 $
% $Revision: 1.11 $
%
% $Log: mcf.m,v $
% Revision 1.11  2009/03/24 08:32:57  marwan
% copyright address changed
%
% Revision 1.10  2007/12/20 16:26:06  marwan
% changed gpl splash behaviour
%
% Revision 1.9  2005/03/16 11:19:02  marwan
% help text modified
%
% Revision 1.8  2004/11/10 07:05:43  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

error(nargchk(2,4,nargin))
error(nargoutchk(0,2,nargout))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% splash the GPL

filename='ace';
which_res=which([filename,'.m']);
gplrc_path=[strrep(which_res,[filename,'.m'],''), 'private'];
gplrc_file=[gplrc_path, filesep, '.gpl.',filename];
if ~exist(gplrc_path)
  mkdir(strrep(which_res,[filename,'.m'],''),'private')
end
if ~exist(gplrc_file)
  fid=fopen(gplrc_file,'w');
  fprintf(fid,'%s\n','If you delete this file, the GNU Public License will');
  fprintf(fid,'%s','splash up at the next time the programme starts.');
  fclose(fid);

  if exist('gpl')
    txt=gpl;
  else
    txt={'The GNU General Public License was not found.'};
  end
  h=figure('NumberTitle','off',...,
         'ButtonDownFcn','close',...
         'Name','GNU General Public License');
  ha=get(h,'Position');
  h=uicontrol('Style','Listbox',...
            'ButtonDownFcn','close',...
            'CallBack','close',...
            'Position',[0 0 ha(3) ha(4)],...
	    'FontName','Courier',...
	    'BackgroundColor',[.8 .8 .8],...
	    'String',txt);
  waitfor(h)
end

if size(y,1)<size(y,2)
  y=y';
end
if size(x,1)<size(x,2)
  x=x';
end

if nargin < =3
  T=round(.02*size(y,1));
end
if nargin < =2  | isempty(w)
  w=round(.1*size(y,1));w=5;
end

if T>0, h=waitbar(0,'Please wait ...','Name','Delay progress','ShowTime',0); end
for i=0:T;
if T>0 & i>0, waitbar(i/(T*2)), end
a=x(i+1:end-T+i);
b=y(1:end-T,:);
c(T-i+1)=ace(a,b,w);
end

for i=1:T;
waitbar((T+i)/(T*2))
a=x(1:end-T);
b=y(0+i:end-T+i-1,:);
c(i+T+1)=ace(a,b,w);
end
if T>0, close(h), end

i=-T:T;

if ~nargout
plot(i,c)
xlabel('Lag')
ylabel('Maximal Correlation \Psi')
axis([-T T 0 1])
line([0 0], [0 1],'LineStyle',':','Color',[0 0 0])

[c1 c2]=max(c);
tx=['Maximum: ', num2str(c1), ' at ' ,num2str(i(c2))];
title(tx)
end

if nargout==1
  a_out=c;
elseif nargout==2
  a_out=i;
  b_out=c;
end
