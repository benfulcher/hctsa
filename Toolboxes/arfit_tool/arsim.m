function [v]=arsim(w,A,C,n_ntr,ndisc)
%ARSIM	Simulation of AR process.	
%
%  v=ARSIM(w,A,C,n) simulates n time steps of the AR(p) process
%
%     v(k,:)' = w' + A1*v(k-1,:)' +...+ Ap*v(k-p,:)' + eta(k,:)', 
%
%  where A=[A1 ... Ap] is the coefficient matrix, and w is a vector of
%  intercept terms that is included to allow for a nonzero mean of the
%  process. The vectors eta(k,:) are independent Gaussian noise
%  vectors with mean zero and covariance matrix C.
%
%  The p vectors of initial values for the simulation are taken to
%  be equal to the mean value of the process. (The process mean is
%  calculated from the parameters A and w.) To avoid spin-up effects,
%  the first 10^3 time steps are discarded. Alternatively,
%  ARSIM(w,A,C,n,ndisc) discards the first ndisc time steps.
%
%  ARSIM(w,A,C,[n, ntr]) generates ntr realizations (trials) of
%  length n of the AR(p) process, which are output as the matrices
%  v(:,:,itr) with itr=1,...,ntr. 
  
%  Modified 13-Oct-00
%  Author: Tapio Schneider
%          tapio@gps.caltech.edu

  m       = size(C,1);                  % dimension of state vectors 
  p       = size(A,2)/m;                % order of process
  n       = n_ntr(1);                   % number of time steps
  if size(n_ntr) == 1
    ntr   = 1;
  else
    ntr   = n_ntr(2);
  end
  
  if (p ~= round(p)) 
    error('Bad arguments.'); 
  end

  if (length(w) ~= m | min(size(w)) ~= 1)
    error('Dimensions of arguments are mutually incompatible.')
  end 
  w       = w(:)';                      % force w to be row vector

  % Check whether specified model is stable
  A1 	  = [A; eye((p-1)*m) zeros((p-1)*m,m)];
  lambda  = eig(A1);
  if any(abs(lambda) > 1)
    warning('The specified AR model is unstable.')
  end
  
  % Discard the first ndisc time steps; if ndisc is not given as input
  % argument, use default
  if (nargin < 5) 
    ndisc = 10^3; 
  end
  
  % Compute Cholesky factor of covariance matrix C
  [R, err]= chol(C);                    % R is upper triangular
  if err ~= 0
    error('Covariance matrix not positive definite.')
  end
    
  % Get ntr realizations of ndisc+n independent Gaussian
  % pseudo-random vectors with covariance matrix C=R'*R
  for itr=1:ntr
    randvec(:, :, itr) = randn([ndisc+n,m])*R;
  end

  % Add intercept vector to random vectors
  randvec = randvec + repmat(w, [ndisc+n, 1, ntr]);
  
  % Get transpose of system matrix A (use transpose in simulation because 
  % we want to obtain the states as row vectors)
  AT      = A';

  % Take the p initial values of the simulation to equal the process mean, 
  % which is calculated from the parameters A and w
  if any(w)
    %  Process has nonzero mean    mval = inv(B)*w'    where 
    %             B = eye(m) - A1 -... - Ap; 
    %  Assemble B
    B 	 = eye(m);
    for j=1:p
      B = B - A(:, (j-1)*m+1:j*m);
    end
    %  Get mean value of process
    mval = w / B';

    %  The optimal forecast of the next state given the p previous
    %  states is stored in the vector x. The vector x is initialized
    %  with the process mean.
    x    = repmat(mval, [p, 1]);
  else
    %  Process has zero mean
    x    = zeros(p,m); 
  end
  
  % Initialize state vectors
  u      = repmat([x; zeros(ndisc+n,m)], [1, 1, ntr]);
  
  % Simulate ntr realizations of n+ndisc time steps. In order to be
  % able to make use of Matlab's vectorization capabilities, the
  % cases p=1 and p>1 must be treated separately.
  if p==1
    for itr=1:ntr
      for k=2:ndisc+n+1; 
        x(1,:) = u(k-1,:,itr)*AT;
        u(k,:,itr) = x + randvec(k-1,:,itr);
      end
    end
  else
    for itr=1:ntr
      for k=p+1:ndisc+n+p; 
        for j=1:p;
          x(j,:) = u(k-j,:,itr)*AT((j-1)*m+1:j*m,:);
        end
        u(k,:,itr) = sum(x)+randvec(k-p,:,itr);
      end
    end
  end
  
  % return only the last n simulated state vectors
  v = u(ndisc+p+1:ndisc+n+p,:,:); 





