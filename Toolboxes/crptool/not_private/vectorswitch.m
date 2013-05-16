function vectorswitch
% VECTORSWITCH   Switches the sign of phase space vectors
%    Used by CRP Toolbox

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
% $Log: vectorswitch.m,v $
% Revision 4.8  2009/03/24 08:36:04  marwan
% copyright address updated
%
% Revision 4.7  2004/11/10 07:04:30  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

global props

	ds = get(findobj('Tag','Dim','Parent',gcf),'UserData');
	j=str2num(get(gco,'String'));
	if ds(j,j)==1
	   set(gco,props.buttonActive,'FontWeight','bold');
	   ds(j,j)= -1;
	else
	   set(gco,props.button,'FontWeight','normal');
	   ds(j,j)= 1;
	end
	h=findobj('Tag','Dim');
	set(h(1),'UserData', ds);
