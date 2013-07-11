function x=mov_av(y,k);

% function x=mov_av(y,k);
%
% subtracts, from a time series y a moving avegrage consisting of the
% last k terms.
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk

if nargin < 2
	k=5;
end;

for i=k+1:length(y)
	x(i-k)=y(i)-mean(y(i-k:i-1));
end;
