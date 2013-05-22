function out = EN_sampenc2(y,M,r,preprocess)

% Ben Fulcher 11/09
% Added out statistics and optional preprocessing input

if nargin > 3 % specified a preprocessing for y
    if strcmp(preprocess,'diff1')
        % this produces 'control entropy'
        y = zscore(diff(y));
    end
end

% function [e,A,B]=sampenc(y,M,r);
%
% Input
% 
% y input data
% M maximum template length
% r matching tolerance
% 
% Output
% 
% e sample entropy estimates for m=0,1,...,M-1
% A number of matches for m=1,...,M
% B number of matches for m=1,...,M excluding last point

N = length(y);
lastrun = zeros(1,N);
run = zeros(1,N);
A = zeros(M,1);
B = zeros(M,1);
p = zeros(M,1);
e = zeros(M,1);
for i = 1:(N-1)
   nj = N-i;
   y1 = y(i);
   for jj = 1:nj
      j = jj + i;      
      if abs(y(j)-y1)<r
         run(jj) = lastrun(jj)+1;
         M1 = min(M,run(jj));
         for m = 1:M1           
            A(m) = A(m)+1;
            if j < N
               B(m) = B(m)+1;
            end
         end
      else
         run(jj) = 0;
      end
   end
   for j = 1:nj
      lastrun(j) = run(j);
   end
end

NN = N*(N-1)/2;
p(1) = A(1)/NN;
e(1) = -log(p(1));

for m = 2:M
   p(m) = A(m)/B(m-1);
   e(m) = -log(p(m));
end


%% Give outputs
for i = 1:M
    eval(['out.p' num2str(i) ' = p(' num2str(i) ');']);
    eval(['out.sampen' num2str(i) ' = e(' num2str(i) ');']);
end

out.meanchsampen = mean(diff(e));
out.meanchp = mean(diff(p));

end