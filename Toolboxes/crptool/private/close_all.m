function close_all
% CLOSE_ALL   Closes all CRP windows.
%    Used by CRP Toolbox

% Copyright (c) 1998-2003 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2004/11/10 07:04:28 $
% $Revision: 4.7 $
%
% $Log: close_all.m,v $
% Revision 4.7  2004/11/10 07:04:28  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

global  props

  if ~isempty(findobj('Tag','CRPFig')), delete(findobj('Tag','CRPFig')), end
  root_ud=get(0,'UserData'); 
  if isstruct(root_ud)
    if isfield(root_ud,'crp')
      root_ud=rmfield(root_ud,'crp');
      if length(fieldnames(root_ud))==1
        if isfield(root_ud,'old'); root_ud=root_ud.old; end
      end
    end
  end
  try, set(0,'UserData',root_ud,props.root), end
  
  clear all
  disp('Thank you for using CRP toolbox.')
