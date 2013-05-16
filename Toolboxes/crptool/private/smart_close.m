function smart_close(hCRP,hCtrl)
% SMART_CLOSE   Closes the current CRP window.
%    Used by CRP Toolbox

% Copyright (c) 1998-2003 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2004/11/10 07:04:29 $
% $Revision: 4.7 $
%
% $Log: smart_close.m,v $
% Revision 4.7  2004/11/10 07:04:29  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

global props

[h h1]=strtok(get(hCRP,'Name'),'(');
if ~isempty(h1)
  h1([1, end])=[];
  if ishandle(hCtrl), delete(hCtrl), end
  if ishandle(hCRP), delete(hCRP), end
  root_ud=get(0,'UserData'); 
  if isstruct(root_ud)
    if isfield(root_ud,'crp')
      root_ud.crp(root_ud.crp==str2num(h1))=[];
      if isempty(root_ud.crp), root_ud=rmfield(root_ud,'crp'); end
    end
    if length(fieldnames(root_ud))==1
      if isfield(root_ud,'old'); root_ud=root_ud.old; end
    end
  end
  try, set(0,'UserData',root_ud,props.root), end
  
  clear all
end
