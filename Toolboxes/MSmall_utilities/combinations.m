function Z=combinations(a,b);

% function Z=combinations(a,b);
% or function Z=combinations(X);
% 
% Produces an a^b by b matrix, each row of which is one of the a^b possible
% elements of the set 
%     { (x1,x2,x3,x3,...,xb) : xi is in {0,1,2,...,(a-1)} for i=1,2,...b }
%
% I think it does it in the quickest possiable way. (The calculation is 
% vectorised).
% 
% If there is one input arguement then [a,b]=size(X).
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk

if nargin < 2
	[a,b]=size(a);
end;

% I'm not even going to try to explain whats going on below: it just works!
i=1:(a^b);
j=1:b;
ii=(i'*ones(size(j)))'; jj=(j'*ones(size(i)))';
Z=(ii-1)-fix((i-1)'*(1./(a.^(j+1))))'.*(a.^(jj+1))';
Z=mod(fix(Z'./(a.^(jj-1))),a);
