function out = rp_plugin
% Cross Recurrence Plot Plugin
%    This binary will be called from some programmes of CRP 
%    toolbox for Matlab.
%
%    Please ensure that you have a compatible version for 
%    your platform/ operating system.

% Copyright (c) 2005 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2005/04/19 $
% $Revision: 1.3 $
%
% $Log: rp_plugin.m,v $
%
%
% This program is part of the new generation XXII series.
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

plugin_path = fileparts(which(mfilename));

if ispc
    plugin_name = 'rp_plugin.exe';
else
    plugin_name = 'rp_plugin';
    if isunix, fileattrib([plugin_path,filesep,plugin_name], '+x'), end
end

if nargout == 1
    out = plugin_name;
else
    try
          [status plugin_version] = system([plugin_path,filesep,plugin_name,' -V']);
    catch
          status = 1;
    end
    
    crp_installed = help('CRPtool');
    if isempty(crp_installed)
        disp('CRP Toolbox not installed! This plugin works only within the')
        disp('toolbox environment. Please install the CRP Toolbox first.')
        return
    end
    
    
    if status
        disp('This plugin doesn''t work for your system. Please ensure that you')
        disp('have a compatible version for your platform/ operating system.')
    else
        disp(plugin_version)
        disp('This is a plugin for the CRP toolbox. It is called automatically')
        disp('from some programmes of this toolbox.')
    end
    
end

