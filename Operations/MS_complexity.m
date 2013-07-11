function out = MS_complexity(x,n,preproc)
% Wrapper code around Michael Small's original 'complexity.m' code. Uses
% the mex file complexitybs.
% http://small.eie.polyu.edu.hk/matlab/
% Ben Fulcher 19/2/2010

if nargin < 2 || isempty(n)
    n=2;
end
if nargin < 3
    preproc = [];
end

if strcmp(preproc,'diff')
    x = benzscore(diff(x));
end

%_______________________________________________________________
% cmp = complexity(x,n);
%
% calculate the Lempel-Ziv complexity of the n-bit encoding of x. 
%
% cmp is the normalised complexity, that is the number of distinct
% symbol sequences in x, divided by the expected number of distinct 
% symbols for a noise sequence.
%
% Algorithm is implemented in complexitybs.c
%
% M. Small
% ensmall@polyu.edu.hk
% 7/10/04

if length(n)>1
   for ni = 1:length(n),
        cmp(ni) = complexity(x,n(ni));
   end
else
    %do the binning, with equi-probably bins
    x=x(:);
    nx=length(x);
    [xn,xi]=sort(x);
%     y=zeros(nx,1);
    y = 1:nx;
    y = floor(y.*(n/(nx+1)));
    x(xi) = y;

    %compute complexity with complexitybs
    cmp = complexitybs(x+eps);
end

%______________________________________________________________
% Ok, so we have cmp, the complexity
out = cmp;

end