function change_colormap(hCtrl,hCRP,cm)
% CHANGE_COLORMAP   Changes the current colormap.

% Copyright (c) 1998-2003 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2004/11/10 07:04:28 $
% $Revision: 4.7 $
%
% $Log: change_colormap.m,v $
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
