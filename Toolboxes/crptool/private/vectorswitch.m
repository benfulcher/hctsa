function vectorswitch
% VECTORSWITCH   Switches the sign of phase space vectors
%    Used by CRP Toolbox

% Copyright (c) 1998-2003 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2004/11/10 07:04:30 $
% $Revision: 4.7 $
%
% $Log: vectorswitch.m,v $
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
