function out = EN_sampenc(y,M,r)
%function [e,A,B]=sampenc(y,M,r);
%
%Input
%
%y input data
%M maximum template length
%r matching tolerance
%
%Output
%
%e sample entropy estimates for m=0,1,...,M-1
%A number of matches for m=1,...,M
%B number of matches for m=1,...,M excluding last point

n=length(y);
lastrun=zeros(1,n);
run=zeros(1,n);
A=zeros(M,1);
B=zeros(M,1);
% p=zeros(M,1);
% e=zeros(M,1);
for i=1:(n-1)
   nj=n-i;
   y1=y(i);
   for jj=1:nj
      j=jj+i;      
      if abs(y(j)-y1)<r
         run(jj)=lastrun(jj)+1;
         M1=min(M,run(jj));
         for m=1:M1           
            A(m)=A(m)+1;
            if j<n
               B(m)=B(m)+1;
            end            
         end
      else
         run(jj)=0;
      end      
   end
   for j=1:nj
      lastrun(j)=run(j);
   end
end
if M == 1;
    N=n*(n-1)/2;
    p(1)=A(1)/N;
    e(1)=-log(p(1));
else
    % for m=2:M
    p=A(M)/B(M-1);
    e(M)=-log(p);
    % end
end
out=e(M);