function [X, matext] = crp_plugin(x, y, m, t, e, mflag, hCRP, plugin_path, silent)
% CRP_PLUGIN   Loads and executes the extern plugin
%    Used by CRP Toolbox

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 2005-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:36:04 $
% $Revision: 4.9 $
%
% $Log: crp_plugin.m,v $
% Revision 4.9  2009/03/24 08:36:04  marwan
% copyright address updated
%
% Revision 4.8  2007/05/22 13:49:45  marwan
% added fixed RR support
%
% Revision 4.7  2007/05/22 13:43:45  marwan
% added fixed RR support
%
% Revision 4.6  2006/10/24 14:16:26  marwan
% minor change: sigma in title line of RP shown only for normalised data
%
% Revision 4.5  2006/03/29 13:07:40  marwan
% problems regarding OPRPs and embedding resolved
%
% Revision 4.4  2006/02/14 11:46:53  marwan
% dimension and delay for plugin released
%
% Revision 4.3  2006/02/06 14:44:36  marwan
% plugin support for order patterns
%
% Revision 4.2  2005/04/15 09:03:03  marwan
% minor bugfix in plugin section
%
% Revision 4.1  2005/04/08 09:03:53  marwan
% plugin added
%
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

global errcode nonorm

X = []; matext = '';

try

      while 1
          tmp_xdatafile = tempname;
          if ~exist(tmp_xdatafile), break, end
      end
      while 1
          tmp_ydatafile = tempname;
          if ~exist(tmp_ydatafile), break, end
      end
      while 1
          tmp_rpdatafile = tempname;
          if ~exist(tmp_rpdatafile), break, end
      end
          
      if ~silent, set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Export Data'),drawnow, end
      save(tmp_xdatafile,'x','-ascii','-tabs');
      if ~isequal(x,y), save(tmp_ydatafile,'y','-ascii','-tabs'); end
      
      % prepare the arguments for the external plugin
      if ~silent, set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Compute Recurrence Points'),drawnow, end
      m_str = {'MAX', 'EUC', 'MIN', 'NR', 'RR', 'FAN', 'IN', 'OM', 'OP', 'EUC'};
      dis_sign = 1; if mflag == 10, dis_sign = -1; end

      unix_str = [plugin_path,filesep,rp_plugin,' -m ',num2str(m), ...
                                    ' -t ',num2str(t), ...
                                    ' -e ',num2str(dis_sign * e), ...
                                    ' -n ',m_str{mflag}, ...
                                    ' -w ',num2str(0), ...
                                    ' -i ',tmp_xdatafile];
      if ~isequal(x,y)
          unix_str = [unix_str,' -j ',tmp_ydatafile];
      end
      
      unix_str = [unix_str,' -r ',tmp_rpdatafile, ...
                                    ' -f TIF', ...
                                    ' -s'];
      
      % call the external plugin
      system(unix_str);
      
      % import rp
      if ~silent, set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Import Recurrence Points'),drawnow, end
      X = imread(tmp_rpdatafile);
      X = (double(X)/double(max(X(:))));

      delete(tmp_rpdatafile);
      delete(tmp_xdatafile);
      if ~isequal(x,y), delete(tmp_ydatafile); end
      
      if ~silent, set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'String','Plot Recurrence Points'),drawnow, end

      if ~nonorm, unit = ''; else unit = '\sigma'; end

      switch mflag
      %%%%%%%%%%%%%%%%% maximum norm
        case 1
          errcode=111;
          matext=[num2str(round(100*e)/100) unit ' (fixed distance maximum norm)'];
      %%%%%%%%%%%%%%%%% euclidean norm
        case 2
          errcode=112;
          matext=[num2str(round(100*e)/100) unit ' (fixed distance euclidean norm)'];
      %%%%%%%%%%%%%%%%% minimum norm
        case 3
          errcode=113;
          matext=[num2str(round(100*e)/100) unit ' (fixed distance minimum norm)'];
      %%%%%%%%%%%%%%%%% order patterns
        case 9
          errcode=118;
          matext='';
      %%%%%%%%%%%%%%%%% distance matrix
        case 10
          errcode=119;
          X = X * ( max([max(x), max(y)]) - min([min(x), min(y)]) );
          matext='';
      end

catch   
  warning off
  delete(tmp_rpdatafile);
  delete(tmp_xdatafile);
  delete(tmp_ydatafile);
  warning on
end
