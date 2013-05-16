% Performs discrete total variation denoising (TVD) using a primal-dual
% interior-point solver. It minimizes the following discrete functional:
%
%  E=(1/2)||y-x||_2^2+lambda*||Dx||_1,
%
% over the variable x, given the input signal y, according to each
% value of the regularization parameter lambda > 0. D is the first
% difference matrix. Uses hot-restarts from each value of lambda to speed
% up convergence for subsequent values: best use of this feature is made by
% ensuring that the chosen lambda values are close to each other.
%
% Usage:
% [x, E, s] = l1pwc(y, lambda, display, stoptol, maxiter)
%
% Input arguments:
% - y          Original signal to denoise, size N x 1.
% - lambda     A vector of positive regularization parameters, size L x 1.
%              TVD will be applied to each value in the vector.
% - display    (Optional) Set to 0 to turn off progress display, 1 to turn
%              on. If not specifed, defaults to progress display on.
% - stoptol    (Optional) Precision as determined by duality gap tolerance,
%              if not specified, defaults to 1e-3.
% - maxiter    (Optional) Maximum interior-point iterations, if not
%              specified defaults to 60.
%
% Output arguments:
% - x          Denoised output signal for each value of lambda, size N x L.
% - E          Objective functional at minimum for each lambda, size L x 1.
% - s          Optimization result, 1 = solved, 0 = maximum iterations
%              exceeded before reaching duality gap tolerance, size L x 1.
% - lambdamax  Maximum value of lambda for the given y. If
%              lambda >= lambdamax, the output is the trivial constant
%              solution x = mean(y).
%
% (c) Max Little, 2010. Based around code originally written by 
% S.J. Kim, K. Koh, S. Boyd and D. Gorinevsky. If you use this code for
% your research, please cite:
% M.A. Little, Nick S. Jones (2010)
% "Sparse Bayesian Step-Filtering for High-Throughput Analysis of Molecular
% Machine Dynamics", in 2010 IEEE International Conference on Acoustics,
% Speech and Signal Processing, 2010, ICASSP 2010 Proceedings.
%
% This code is released under the terms of GNU General Public License as
% published by the Free Software Foundation; version 2 or later.

function [x, E, s, lambdamax] = l1pwc(y, lambda, display, stoptol, maxiter)

error(nargchk(2,5,nargin));
if (nargin < 3)
    display = 1;
end
if (nargin < 4)
    stoptol = 1e-3;
end
if (nargin < 5)
    maxiter = 60;
end

y = y(:);

% Search tuning parameters
ALPHA     = 0.01;   % Backtracking linesearch parameter (0,0.5]
BETA      = 0.5;    % Backtracking linesearch parameter (0,1)
MAXLSITER = 20;     % Max iterations of backtracking linesearch
MU        = 2;      % t update

N = length(y);    % Length of input signal y
M = N-1;          % Size of Dx

% Construct sparse operator matrices
I1 = speye(M,M);
O1 = spalloc(M,1,M);
D = [I1 O1]-[O1 I1];

DDT = D*D';
Dy  = D*y;

% Find max value of lambda
lambdamax = max(abs(DDT\Dy));

if (display)
    fprintf('lambda_max=%5.2e\n', lambdamax);
end

L = length(lambda);
x = zeros(N, L);
s = zeros(L, 1);
E = zeros(L, 1);

% Optimization variables set up once at the start
z    = zeros(M,1);   % Dual variable
mu1  = ones(M,1);    % Dual of dual variable
mu2  = ones(M,1);    % Dual of dual variable

% Work through each value of lambda, with hot-restart on optimization
% variables
for l = 1:L
    
    t    =  1e-10; 
    step =  Inf;
    f1   =  z-lambda(l);
    f2   = -z-lambda(l);

    % Main optimization loop
    s(l) = 1;

    if (display)
        fprintf('Solving for lambda=%5.2e, lambda/lambda_max=%5.2e\nIter# Primal    Dual      Gap\n', ...
            lambda(l), lambda(l)/lambdamax);
    end
    for iters = 0:maxiter

        DTz  = (z'*D)';
        DDTz = D*DTz;
        w    = Dy-(mu1-mu2);

        % Calculate objectives and primal-dual gap
        pobj1 = 0.5*w'*(DDT\w)+lambda(l)*sum(mu1+mu2);
        pobj2 = 0.5*DTz'*DTz+lambda(l)*sum(abs(Dy-DDTz));
        pobj = min(pobj1,pobj2);
        dobj = -0.5*DTz'*DTz+Dy'*z;
        gap  = pobj - dobj;

        if (display)
            fprintf('%5d %7.2e %7.2e %7.2e\n', iters, pobj, dobj, gap);
        end

        % Test duality gap stopping criterion
        if (gap <= stoptol)
            s(l) = 1;
            break;
        end

        if (step >= 0.2)
            t = max(2*M*MU/gap, 1.2*t);
        end

        % Do Newton step
        rz      =  DDTz - w;
        S       =  DDT-sparse(1:M,1:M,mu1./f1+mu2./f2);
        r       = -DDTz + Dy + (1/t)./f1 - (1/t)./f2;
        dz      =  S\r;
        dmu1    = -(mu1+((1/t)+dz.*mu1)./f1);
        dmu2    = -(mu2+((1/t)-dz.*mu2)./f2);

        resDual = rz;
        resCent = [-mu1.*f1-1/t; -mu2.*f2-1/t];
        residual= [resDual; resCent];

        % Perform backtracking linesearch
        negIdx1 = (dmu1 < 0); 
        negIdx2 = (dmu2 < 0);
        step = 1;
        if (any(negIdx1))
            step = min( step, 0.99*min(-mu1(negIdx1)./dmu1(negIdx1)) );
        end
        if (any(negIdx2))
            step = min( step, 0.99*min(-mu2(negIdx2)./dmu2(negIdx2)) );
        end

        for liter = 1:MAXLSITER
            newz    =  z  + step*dz;
            newmu1  =  mu1 + step*dmu1;
            newmu2  =  mu2 + step*dmu2;
            newf1   =  newz - lambda(l);
            newf2   = -newz - lambda(l);

            % Update residuals
            newResDual  = DDT*newz - Dy + newmu1 - newmu2;
            newResCent  = [-newmu1.*newf1-1/t; -newmu2.*newf2-1/t];
            newResidual = [newResDual; newResCent];

            if ((max(max(newf1),max(newf2)) < 0) && ...
                    (norm(newResidual) <= (1-ALPHA*step)*norm(residual)))
                break;
            end
            step = BETA*step;
        end

        % Update primal and dual optimization parameters
        z = newz;
        mu1 = newmu1;
        mu2 = newmu2;
        f1 = newf1;
        f2 = newf2;
    end

    x(:,l) = y-D'*z;
    E(l) = 0.5*sum((y-x(:,l)).^2)+lambda(l)*sum(abs(D*x(:,l)));

    % We may have a close solution that does not satisfy the duality gap
    if (iters >= maxiter)
        s(l) = 0;
    end
    
    if (display)
        if (s(l))
            fprintf('Solved to precision of duality gap %5.2e\n', gap);
        else
            fprintf('Max iterations exceeded - solution may be inaccurate\n');
        end
    end

end
