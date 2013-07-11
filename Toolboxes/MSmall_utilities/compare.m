function rms=compare(y,n,m);

%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk

if nargin < 3
	if nargin < 2
		n=1; m=60;
	else
		m=n; n=1;
	end;
end;

y=y-mean(y);
y=y/std(y);
y=y+10;

[a,b]=size(y);
if b>a, y=y'; a=b; end;

for i=n:m
%for j=(i+1):m
	%X=[y((j+1-i):(a-i)) y(1:(a-j))];
	%z=y((j+1):a);
	X=y(1:(a-i)); z=y((i+1):a);
	lam=X\z;
	e=X*lam-z;
	%rms(i,j)=RMS(e);
	rms(i)=RMS(e);
%end; disp(int2str(i));
end;

plot(rms); grid

