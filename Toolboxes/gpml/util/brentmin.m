%% BRENTMIN: Brent's minimization method in one dimension
function [xmin,fmin,funccount,varargout] = ...
                                  brentmin(xlow,xupp,Nitmax,tol,f,nout,varargin)
% code taken from
%    § 10.2 Parabolic Interpolation and Brent's Method in One Dimension
%    Press, Teukolsky, Vetterling & Flannery
%    Numerical Recipes in C, Cambridge University Press, 2002
%
% [xmin,fmin,funccout,varargout] = BRENTMIN(xlow,xupp,Nit,tol,f,nout,varargin)
%    Given a function f, and given a search interval this routine isolates 
%    the minimum of fractional precision of about tol using Brent's method.
% 
% INPUT
% -----
% xlow,xupp:  search interval such that xlow<=xmin<=xupp
% Nitmax:     maximum number of function evaluations made by the routine
% tol:        fractional precision 
% f:          [y,varargout{:}] = f(x,varargin{:}) is the function
% nout:       no. of outputs of f (in varargout) in addition to the y value
%
% OUTPUT
% ------
% fmin:      minimal function value
% xmin:      corresponding abscissa-value
% funccount: number of function evaluations made
% varargout: additional outputs of f at optimum
%
% Copyright (c) by Hannes Nickisch 2010-01-10.

if nargin<6, nout = 0; end
varargout = cell(nout,1);

% tolerance is no smaller than machine's floating point precision
tol = max(tol,eps);

% Evaluate endpoints
fa = f(xlow,varargin{:});
fb = f(xupp,varargin{:});
funccount = 2; % number of function evaluations
% Compute the start point
seps = sqrt(eps);
c = 0.5*(3.0 - sqrt(5.0));% golden ratio
a = xlow; b = xupp;
v = a + c*(b-a);
w = v; xf = v;
d = 0.0; e = 0.0;
x = xf; [fx,varargout{:}] = f(x,varargin{:});
funccount = funccount + 1;

fv = fx; fw = fx;
xm = 0.5*(a+b);
tol1 = seps*abs(xf) + tol/3.0;
tol2 = 2.0*tol1;

% Main loop
while ( abs(xf-xm) > (tol2 - 0.5*(b-a)) )
    gs = 1;
    % Is a parabolic fit possible
    if abs(e) > tol1
        % Yes, so fit parabola
        gs = 0;
        r = (xf-w)*(fx-fv);
        q = (xf-v)*(fx-fw);
        p = (xf-v)*q-(xf-w)*r;
        q = 2.0*(q-r);
        if q > 0.0,  p = -p; end
        q = abs(q);
        r = e;  e = d;

        % Is the parabola acceptable
        if ( (abs(p)<abs(0.5*q*r)) && (p>q*(a-xf)) && (p<q*(b-xf)) )

            % Yes, parabolic interpolation step
            d = p/q;
            x = xf+d;

            % f must not be evaluated too close to ax or bx
            if ((x-a) < tol2) || ((b-x) < tol2)
                si = sign(xm-xf) + ((xm-xf) == 0);
                d = tol1*si;
            end
        else
            % Not acceptable, must do a golden section step
            gs=1;
        end
    end
    if gs
        % A golden-section step is required
        if xf >= xm, e = a-xf;    else e = b-xf;  end
        d = c*e;
    end

    % The function must not be evaluated too close to xf
    si = sign(d) + (d == 0);
    x = xf + si * max( abs(d), tol1 );
    [fu,varargout{:}] = f(x,varargin{:});
    funccount = funccount + 1;

    % Update a, b, v, w, x, xm, tol1, tol2
    if fu <= fx
        if x >= xf, a = xf; else b = xf; end
        v = w; fv = fw;
        w = xf; fw = fx;
        xf = x; fx = fu;
    else % fu > fx
        if x < xf, a = x; else b = x; end
        if ( (fu <= fw) || (w == xf) )
            v = w; fv = fw;
            w = x; fw = fu;
        elseif ( (fu <= fv) || (v == xf) || (v == w) )
            v = x; fv = fu;
        end
    end
    xm = 0.5*(a+b);
    tol1 = seps*abs(xf) + tol/3.0; tol2 = 2.0*tol1;

    if funccount >= Nitmax        
        % typically we should not get here
        % warning(sprintf(['Maximum number of iterations (%d) exceeded:', ...
        %                  'precision is not guaranteed'],Nitmax))
        % fprintf('[%1.3f,%1.3f,%1.3f]\n',xlow,xf,xupp)
        break
    end
end % while

% check that endpoints are less than the minimum found
if ( (fa < fx) && (fa <= fb) )
    xf = xlow; fx = fa;
elseif fb < fx
    xf = xupp; fx = fb;
end
fmin = fx;
xmin = xf;
