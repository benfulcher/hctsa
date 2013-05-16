function switch_unthresholded(hCRP)
% SWITCH_UNTHRESHOLDED   Switches between thesholded/unthesholded CRP plot
%    Used by CRP Toolbox

% Copyright (c) 1998-2004 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2008/02/25 11:45:27 $
% $Revision: 4.9 $
%
% $Log: switch_unthresholded.m,v $
% Revision 4.9  2008/02/25 11:45:27  marwan
% fix of the colorbar bug
%
% Revision 4.8  2005/12/09 19:19:17  marwan
% && and || changed to & and | (upwards compatibility)
%
% Revision 4.7  2004/11/10 07:04:29  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.


  if get(gcbo,'value')==12 & strcmp(get(gcbo,'tag'),'Method')
  set(findobj('Tag','Unthresh','Parent',gcf),'enable','off')
  set(findobj('Tag','Dim','Parent',gcf),'enable','off')
  set(findobj('Tag','Dimtext','Parent',gcf),'enable','off')
  set(findobj('Tag','Delay','Parent',gcf),'enable','off')
  set(findobj('Tag','Delaytext','Parent',gcf),'enable','off')
  return
  end
  
  if get(gcbo,'value')~=12 & strcmp(get(gcbo,'tag'),'Method')
  set(findobj('Tag','Unthresh','Parent',gcf),'enable','on')
  set(findobj('Tag','Dim','Parent',gcf),'enable','on')
  set(findobj('Tag','Dimtext','Parent',gcf),'enable','on')
  set(findobj('Tag','Delay','Parent',gcf),'enable','on')
  set(findobj('Tag','Delaytext','Parent',gcf),'enable','on')
  return
  end

  if get(findobj('Tag','Unthresh','Parent',gcf),'value')==1
  set(findobj('Tag','Size','Parent',gcf),'enable','off')
  set(findobj('Tag','Sizetext','Parent',gcf),'enable','off')
  set(findobj('Tag','Method','Parent',gcf),'enable','off')
  set(findobj('Tag','Log','Parent',gcf),'enable','on')
  set(findobj('Tag','Colorbar','Parent',hCRP),'visible','on')
  set(get(findobj('Tag','Colorbar','Parent',hCRP),'children'),'visible','on')
  end
  if get(findobj('Tag','Unthresh','Parent',gcf),'value')==0
  set(findobj('Tag','Size','Parent',gcf),'enable','on')
  set(findobj('Tag','Sizetext','Parent',gcf),'enable','on')
  set(findobj('Tag','Method','Parent',gcf),'enable','on')
  set(findobj('Tag','Colorbar','Parent',hCRP),'visible','off')
  set(get(findobj('Tag','Colorbar','Parent',hCRP),'children'),'visible','off')
  set(findobj('Tag','Log','Parent',gcf),'enable','off')
  end
