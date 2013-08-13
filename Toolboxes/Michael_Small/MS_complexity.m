% cmp = MS_complexity(x,n);
%
% calculate the Lempel-Ziv complexity of the n-bit encoding of x. 
%
% cmp is the normalised complexity, that is the number of distinct
% symbol sequences in x, divided by the expected number of distinct 
% symbols for a noise sequence.
%
% Algorithm is implemented in complexitybs.c
%
% Michael Small
% michael.small@uwa.edu.au, http://school.maths.uwa.edu.au/~small/
% 7/10/04
% For further details, please see M. Small. Applied Nonlinear Time Series
% Analysis: Applications in Physics, Physiology and Finance. Nonlinear Science
% Series A, vol. 52. World Scientific, 2005. (ISBN 981-256-117-X) and the
% references therein.
% (minor cosmetic changes by Ben Fulcher, 2010)

function cmp = MS_complexity(x,n);

if nargin<2,
    n=2;
end;

if length(n)>1,
   for ni=1:length(n),
        cmp(ni) = MS_complexity(x,n(ni)); 
   end;    
else,
    
    if 1,
        %do the binning, with equi-probably bins
        x=x(:);
        nx=length(x);
        [xn,xi]=sort(x+eps*randn(size(x))); %introduce randomness for ties
        y=zeros(nx,1);
        y=1:nx;
        y=floor(y.*(n/(nx+1)));
        x(xi)=y;
    else,
        %do binning with equal width bins
        x=x(:);
        nx=length(x);
        minx=min(x);
        maxx=max(x);
        stepx=(maxx-minx)/n;
        y=zeros(nx,1);
        while minx<maxx,
            minx=minx+stepx;
            y=y+double(x<minx);
        end;
        x=floor(y);
    end;
    
    %compute complexity with complexitybs
    cmp=MS_complexitybs(x);

end;