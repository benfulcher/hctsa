function makescript(d, filename)

%tstoolbox/@description/makescript
%   Syntax:
%     * makescript (signal, scriptfilename)
%
%   creates a Matlab m-file that contains exactly the the processing steps
%   that have been applied to get the input signal. This gives a kind of
%   macro facility for tstool.
%   Example
%          signal s was calculated through several processing steps from
%          signal s0 (the raw or original signal) Now makescript(s,
%          'foo.m') will create a Matlab m-file named foo.m which, applied
%          to s0, will give s.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,2);

[path,name,ext] = fileparts(filename); 

n = findstr(name, '.m');
if isempty(n)
	funname = name;
else
	funname = name(1:n-1);
end

if length(d.commandlines)==0
   error('Signal has no commands that can be used to create a script');
end

[fid, message] = fopen(filename, 'wt');	% open as text file

if fid == -1
	error(message)
end

line = ['function s = ' funname '(s)\n'];
fprintf(fid, line);

for i=1:length(d.commandlines)
   line = [get(d.commandlines, i) '\n'];
   fprintf(fid, line);
end

fclose(fid);
