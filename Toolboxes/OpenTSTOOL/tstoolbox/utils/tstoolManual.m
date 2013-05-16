function tstoolManual
% Displays HTML OpenTSTOOL manual in help browser

base = fileparts(which('tstoolInit'));
if ~isempty(base)
	web(fullfile(base,'Doc/HTML/index.html'),'-helpbrowser');
else
	error('OpenTSTOOL directories are not in Matlab search path')
end

