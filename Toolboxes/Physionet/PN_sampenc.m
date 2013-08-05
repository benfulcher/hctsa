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

function [e, p, A, B] = PN_sampenc(y,M,r,justM)
% Modified very slightly from original code sampenc.m from 
% http://physionet.org/physiotools/sampen/
% http://www.physionet.org/physiotools/sampen/matlab/1.1/sampenc.m
% Code by DK Lake (dlake@virginia.edu), JR Moorman and Cao Hanqing.
% 
% Added input checking, and altered a few minor things, including
% how the outputs are presented.
% Also added input 'justM', to give e just for the given M, and not for
% all other m up to it
% Ben Fulcher, 2010

% Check inputs
if nargin < 2 || isempty(M)
    M = 1;
end
if nargin < 3 || isempty(r)
    r = 0.1;
end
if nargin < 4 || isempty(justM)
    justM = 0;
end

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

% Calculate for m = 1
NN = N*(N-1)/2;
p(1) = A(1)/NN;
e(1) = -log(p(1));

% Calculate for m > 1, up to M
for m = 2:M
   p(m) = A(m)/B(m-1);
   e(m) = -log(p(m));
end

% output vector p for m = 1, ..., M
% output vector e for m = 1, ..., M

if justM
    % just output the entropy and probability at the maximum m
    e = e(end);
    p = p(end); % (just in case)
end

end