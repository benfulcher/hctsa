function print_error(programme,z,x,y,in,mflag,action)
% PRINT_ERROR   Prints an error message.
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
% $Log: print_error.m,v $
% Revision 4.8  2009/03/24 08:36:04  marwan
% copyright address updated
%
% Revision 4.7  2004/11/10 07:04:29  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

global errcode props

  if ~isempty(findobj('Tag','TMWWaitbar')), delete(findobj('Tag','TMWWaitbar')), end
  if ischar(in), in2=in; else, in2=[]; end
  in=whos('in');
  if ~strcmpi(lasterr,'Interrupt')
    fid=fopen('error.log','w');
    err=fprintf(fid,'%s\n','Please send us the following error report. Provide a brief');
    err=fprintf(fid,'%s\n','description of what you were doing when this problem occurred.');
    err=fprintf(fid,'%s\n','E-mail or FAX this information to us at:');
    err=fprintf(fid,'%s\n','    E-mail:  marwan@pik-potsdam.de');
    err=fprintf(fid,'%s\n','       Fax:  ++49 +331 288 2640');
    err=fprintf(fid,'%s\n\n\n','Thank you for your assistance.');
    err=fprintf(fid,'%s\n',repmat('-',50,1));
    err=fprintf(fid,'%s\n',datestr(now,0));
    err=fprintf(fid,'%s\n',['Matlab ',char(version),' on ',computer]);
    err=fprintf(fid,'%s\n',repmat('-',50,1));
    err=fprintf(fid,'%s\n',x);
    err=fprintf(fid,'%s\n',y);
    if ~isempty(action)
      err=fprintf(fid,'%s\n',[' during ==> ',programme,':',action]);
    else
      err=fprintf(fid,'%s\n',[' during ==> ',programme]);
    end
    if ~isempty(mflag)
      err=fprintf(fid,'%s\n',[' method ==> ',num2str(mflag)]);
    end
    err=fprintf(fid,'%s',[' input ==> ',in.class]);
    if ~isempty(in2), err=fprintf(fid,'\t%s\n',[' (',in2,')']); end
    if ~isempty(errcode)
      err=fprintf(fid,'%s\n',[' errorcode ==> ',num2str(errcode)]);
    else
      err=fprintf(fid,'%s\n',[' errorcode ==> no errorcode available']);
    end
    err=fprintf(fid,'%s\n',' workspace dump ==>');
    if ~isempty(z), 
      err=fprintf(fid,'%s\n',['Name',char(9),'Size',char(9),'Bytes',char(9),'Class']);
    for j=1:length(z);
      err=fprintf(fid,'%s\n',[z(j).name,char(9),num2str(z(j).size),char(9),num2str(z(j).bytes),char(9),z(j).class]);
    end, end
    err=fclose(fid);
    disp('----------------------------');
    disp('       ERROR OCCURED ');
    disp(['   during executing ',programme]);
    disp('----------------------------');
    disp(x);
    if ~isempty(action)
      disp(['   during ',action]);
    end
    if ~isempty(errcode)
      disp(['   errorcode is ',num2str(errcode)]);
    else
      disp('   no errorcode available');
    end
    disp('----------------------------');
    disp('   Please send us the error report. For your convenience, ')
    disp('   this information has been recorded in: ')
    disp(['   ',fullfile(pwd,'error.log')]), disp(' ')
    disp('   Provide a brief description of what you were doing when ')
    disp('   this problem occurred.'), disp(' ')
    disp('   E-mail or FAX this information to us at:')
    disp('       E-mail:  marwan@pik-potsdam.de')
    disp('          Fax:  ++49 +331 288 2640'), disp(' ')
    disp('   Thank you for your assistance.')
    warning('on')
  end
  try, set(0,props.root), end
  set(0,'ShowHidden','Off')
