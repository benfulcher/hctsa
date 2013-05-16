function [a_out, b_out]=pss(varargin)
%PSS   Computes phase space size.
%    [Y Z]=PSS(X [,M,T]) computes the maximal (Y) and the 
%    averaged (Z) phase space diameter of embedded data X,
%    with the embedding parameters dimension M and lag T.
%
%    [Y Z]=PSS(..., NORM) computes the maximal and averaged
%    phase space diameter by using a specified norm, where 
%    NORM can be
%      maxnorm     - Maximum norm.
%      euclidean   - Euclidean norm (default).
%      minnorm     - Minimum norm.
%
%    See also PHASESPACE, CRP, CRQA.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 1998-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:34:24 $
% $Revision: 2.2 $
%
% $Log: pss.m,v $
% Revision 2.2  2009/03/24 08:34:24  marwan
% copyright address changed
%
% Revision 2.1  2004/11/10 07:10:06  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

error(nargchk(1,4,nargin));
if nargout>2, error('Too many output arguments'), end

i_double=find(cellfun('isclass',varargin,'double'));
i_char=find(cellfun('isclass',varargin,'char'));
check_meth={'ma','eu','mi'}; 	% maxnorm, euclidean, minnorm
temp_meth=0;

   if ~isempty(i_char)
      for i=1:length(i_char), 
         varargin{i_char(i)}(3)='0';
         temp_meth=strcmpi(varargin{i_char(i)}(1:2),check_meth');
      end
      method=min(find(temp_meth));
      if isempty(method), method=2; end
      if method>3, method=3; end
   else
      method=2;
   end
   x=double(varargin{i_double(1)});
   if length(i_double)>1 & max(size(varargin{i_double(2)}))==1
     m=varargin{i_double(2)};
   else
     m=1;
   end
   if length(i_double)>2 & max(size(varargin{i_double(3)}))==1
     t=varargin{i_double(3)};
   else
     t=1;
   end

  Nx=length(x);NX=Nx-t*(m-1);
  i=(1:NX)';j=0:t:0+(m-1)*t;
  i=reshape(i(:,ones(1,m))+j(ones(NX,1),:),m*NX,1);
  x1=x(i);
  x2=reshape(x1,NX,m);
  [NX, mx] = size(x2);


  px = permute(x2, [ 1 3 2 ]);
  py = permute(x2, [ 3 1 2 ]);
  s1 = px(:,ones(1,NX),:) - py(ones(1,NX),:,:);

switch(method)
case{1}
  max_dist = max(x(1:end-(m-1)*t))-min(x(1:end-(m-1)*t));
  s = max(abs(s1),[],3);
case{2}
  s = sum(s1.^2, 3);
  max_dist = sqrt(max(s(:)));
case{3}
  s = sum(abs(s1), 3);
  max_dist = max(sum(x2,2))-min(sum(x2,2));
end
mean_dist = sqrt(mean(s(:)));

if nargout==2
  b_out=mean_dist;
end
if nargout>0
  a_out=max_dist;
else
  max_dist
end
