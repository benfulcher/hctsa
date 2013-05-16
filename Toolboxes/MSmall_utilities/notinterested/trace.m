function [level,is_on] = trace(p1,p2)

% [level,is_on] = trace(...)
% set or examine the global trace variable
%
% Setting trace:
%   trace('on') & trace('off')  turn tracing on or off & does not change level
%   trace(level) sets the trace level and does not change whether on or off
%   trace('on',level) etc.. sets level and on/off (args can be any order)
%
% Examining trace
%   0 or 1 outputs return level if tracing is on, [] if off
%   2 outputs returns level and is_on = 1 if tracing is on and 0 if off
%
% Note: When a session starts tracing is off and the level is []
%       This means when tracing is first turned on the level must
%       also be set for tracing to come into effect


% Copyright (c) 1994 by Kevin Judd.  
% Please see the copyright notice included
% in this distribution for full details.
%
% NAME trace.m
%   $Id$

global trace_level;
global trace_is_on;

for i=1:nargin
  p= eval(['p' num2str(i)]);
  if isstr(p)
    trace_is_on= strcmp(p,'on') | strcmp(p,'On');
  else
    trace_level= p;
  end;
end;

if nargout<2
  if trace_is_on
    level= trace_level;
  else
    level=0;is_on=0; %need this for picky matlab5
  end;
else
  level= trace_level;
  is_on= trace_is_on;
end;
