function splash_gpl(filename)
% SPLASH_GPL   Splashs up the Gnu General Public License
%    Used by CRP Toolbox

% Copyright (c) 1998-2003 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2004/11/10 07:04:29 $
% $Revision: 4.7 $
%
% $Log: splash_gpl.m,v $
% Revision 4.7  2004/11/10 07:04:29  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

global props
init_properties

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
         'Name','GNU General Public License',...
	 props.window);
  ha=get(h,'Position');
  h=uicontrol(props.listbox,...
            'ButtonDownFcn','close',...
            'CallBack','close',...
            'Position',[0 0 ha(3) ha(4)],...
	    'String',txt);
  waitfor(h)
end
