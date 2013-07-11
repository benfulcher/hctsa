function [d,i]=nearneigh(X,tau,blocksize)
  
%function [d,i]=nearneigh(X,tau,blocksize)
%
%calculate the nearest (RMS) neighbour of each embedded point
%represented as columns of X.
%tau points either side of each point are excluded (default tau=0);
% i is the index of the nearest neighbours and d are the distances.
%
%nearest neighbours are calculated in a blockwise way in an effort
%to speed up the calculation. This fails horribly, it is quicker to 
%use a simple for loop.
%
%blocksize is the number of points to do at once.
%  
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk
  
if nargin < 2,
  tau=0;
end;

if nargin < 3,
  blocksize=100; %do blocksize points at a time  --- a compromise
		 %between matlabs matrix/memory intensive approach
		 %and CPU ways
end;
	   
[de,nx]=size(X);

%exclude these
exclude=-tau:tau;

b=blocksize+2*tau;
ex=zeros(b);
for i=exclude,
  ex=ex+diag(inf*ones(1,b-abs(i)),i);
end;
ex=ex((tau+1):(end-tau),:);

n=1;
i=[];d=[];
while n<nx,
  %the points
  them=n:(n+blocksize-1);  %do some, but...
  them=them(them<=nx);     %don't go to far
  lt=length(them);
  %and their neighbours
  theirs=(n-tau):(n+blocksize-1+tau);
  these=find(theirs>0 & theirs<=nx);
  theirs=theirs(these);
  
  y=X(:,them);
  y=repmat(y,[1 1 nx]); %create nx copies of each point
  z=repmat(X,[1 1 lt]); %create blocksize copies of X
  z=permute(z,[1 3 2]); %permute z so it matches y
  
  
  %calculate the rms distances
  r=rms(y-z);
  clear y z 
  r=squeeze(r);
  r(1:lt,theirs)=r(1:lt,theirs)+ex(1:lt,these); % to ensure that the
                                           % nearest neighbour is
                                           % on a separate "strand".

  %exclude the end point
  r(1:lt,end)=inf;
  
  %and find the nearest neighbours
  [dst,ind]=min(r');
  ind=ind;%+n-1;
  clear r;

  i=[i ind];
  d=[d dst];
  
  %finally, increment n
  n=n+blocksize;
  
end;
