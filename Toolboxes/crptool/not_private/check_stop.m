function out=check_stop(hCRP,hCtrl,nogui,obj)
% CHECK_STOP   Checks if stop button is pressed.
%    Used by CRP Toolbox

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 1998-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2010/06/29 12:52:08 $
% $Revision: 4.9 $
%
% $Log: check_stop.m,v $
% Revision 4.9  2010/06/29 12:52:08  marwan
% better performance during error debugging
%
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


    global errcode
    errcode=errcode+.001;
    out=0;
    if nogui==0
      h=findobj('Tag','Apply','Parent',hCtrl);
      if strcmpi(get(h(1),'String'),'stopped')
        set(h(1),'String','Apply',...
     	         'ToolTip','Starts the computation - be patient.',...
	         'Callback','crp compute')
        for i=1:length(obj.enable), set(obj.children(i),'Enable',obj.enable{i}); end
        set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')'),'String','Stopped'),drawnow
        setptr([hCRP,hCtrl],'arrow')
        out=1;
      end
    end

