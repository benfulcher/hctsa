%function ind=nearest(x,tau,v);
%
% returns the row vector containing the indicies of the nearest
% neighbours to each of the columns of x. Each point and its tau
% temporal neighbours are excluded from the search. 
% v is an array (not necessarily logical) indicating which columns of x to
% use or, the relative importance of these columns, in the computation
% (i.e. use v(i)*x(i,:), not x(i,:)).  
%
% default : tau=0
%           v=ones(1,length(x(:,1)))
%
% this is a mex version of nearneigh, and provides a speed-up of
% atleast 1000%.
%
% 15/7/04
% Michael Small.
