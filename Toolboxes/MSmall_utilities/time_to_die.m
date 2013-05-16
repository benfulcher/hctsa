function time_to_die(after,before);

% function time_to_die(after);
%
% after is a six element row vector representation of a time.
% if this time is latter than the present time matlab is forced to
% quit.

time=clock;

when=(((after(1)*366+after(2)*31+after(3))*24+after(4))*60+after(5))*60+after(6);
now=(((time(1)*366+time(2)*31+time(3))*24+time(4))*60+time(5))*60+time(6);

if now>when
 disp('Aaahhhhhhhhhh...');
 quit;
end;

