function change_colormap(hCtrl,hCRP,cm)
% CHANGE_COLORMAP   Changes the current colormap.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 1998-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:36:04 $
% $Revision: 4.8 $
%
% $Log: change_colormap.m,v $
% Revision 4.8  2009/03/24 08:36:04  marwan
% copyright address updated
%
% Revision 4.7  2004/11/10 07:04:28  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

  c=str2num(get(gcbo,'Tag'));
  v=ones(1,5);
  if c~=19,
    set(get(get(gcbo,'Parent'),'Children'),'Checked','Off')
    set(gcbo,'Checked','On')
    v=[1 2 4 6 8];
  end
  h2=repmat(cm{c},v(get(findobj('Tag','Log','Parent',hCtrl),'value')),1);
  h1=h2(1:v(get(findobj('Tag','Log','Parent',hCtrl),'value')):end,:);
  set(hCRP,'Colormap',h1)
  
  clear cm cm_old c 
