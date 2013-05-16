function [S, Serr, per, tau, exctn, lambda] = armode(A, C, th)
%ARMODE	Eigendecomposition of AR model.
%
%  [S,Serr,per,tau,exctn]=ARMODE(A,C,th) computes the
%  eigendecomposition of an AR(p) model that has been fitted using
%  ARFIT. The input arguments of ARMODE are output of ARFIT.
%
%  The columns of the output matrix S contain the estimated eigenmodes
%  of the AR model. The output matrix Serr contains margins of error
%  for the components of the estimated eigenmodes S, such that 
%  (S +/- Serr) are approximate 95% confidence intervals for the
%  individual components of the eigenmodes.
%
%  The two-row matrices per and tau contain in their first rows the
%  estimated oscillation period per(1,k) and the estimated damping
%  time tau(1,k) of the eigenmode S(:,k). In their second rows, the
%  matrices per and tau contain margins of error for the periods and
%  damping times, such that 
%     ( per(1,k) +/- per(2,k) )   and   ( tau(1,k) +/- tau(2,k) ) 
%  are approximate 95% confidence intervals for the period and damping
%  time of eigenmode S(:,k).
%  
%  For a purely relaxatory eigenmode, the period is infinite (Inf).
%  For an oscillatory eigenmode, the periods are finite.
%  
%  The excitation of an eigenmode measures its dynamical importance
%  and is returned as a fraction exctn that is normalized such that
%  the sum of the excitations of all eigenmodes equals one.
%
%  See also ARFIT, ARCONF.

%  Modified 13-Oct-00
%  Author: Tapio Schneider
%	   tapio@gps.caltech.edu

  ccoeff   = .95;                       % confidence coefficient
  m 	   = size(C,1);			% dimension of state space
  p 	   = size(A,2) / m; 		% order of model
  if p <= 0 
    error('Order must be greater 0.'); 
  end

  % Assemble coefficient matrix of equivalent AR(1) model
  A1 	   = [A; eye((p-1)*m) zeros((p-1)*m,m)];

  % Eigenvalues and eigenvectors of coefficient matrix of equivalent
  % AR(1) model
  [BigS,d] = eig(A1);  			% columns of BigS are eigenvectors
  lambda   = diag(d);    		% vector containing eigenvalues
  lambda   = lambda(:); 	        % force lambda to be column vector

  % Warning if the estimated model is unstable
  if any(abs(lambda) > 1)
    warning(sprintf(['The estimated AR model is unstable.\n',...
		     '\t Some excitations may be negative.']))
  end
    
  % Fix phase of eigenvectors such that the real part and the
  % imaginary part of each vector are orthogonal
  BigS     = adjph(BigS);

  % Return only last m components of each eigenvector
  S 	   = BigS((p-1)*m+1:p*m, :);

  % Compute inverse of BigS for later use
  BigS_inv = inv(BigS);

  % Recover the matrix Uinv that appears in the asymptotic covariance
  % matrix of the least squares estimator (Uinv is output of AR)
  if (size(th,2) == m*p+1) 
    % The intercept vector has been fitted by AR; in computing
    % confidence intervals for the eigenmodes, this vector is
    % irrelevant. The first row and first column in Uinv,
    % corresponding to elements of the intercept vector, are not
    % needed.
    Uinv   = th(3:size(th,1), 2:size(th,2));

  elseif (size(th,2) == m*p)
    %  No intercept vector has been fitted
    Uinv   = th(2:size(th,1), :);
  else
    error('Input arguments of ARMODE must be output of ARFIT.')
  end
  % Number of degrees of freedom 
  dof 	 = th(1,1);             
  % Quantile of t distribution for given confidence coefficient and dof
  t      = tquant(dof, .5+ccoeff/2); 
  
  % Asymptotic covariance matrix of estimator of coefficient matrix A
  Sigma_A  = kron(Uinv, C);

  % Noise covariance matrix of system of relaxators and oscillators
  CovDcpld = BigS_inv(:, 1:m) * C * BigS_inv(:, 1:m)';

  % For each eigenmode j: compute the period per, the damping time
  % tau, and the excitation exctn; also get the margins of error for
  % per and tau
  for j=1:m*p				% eigenmode number
    a		= real(lambda(j)); 	% real part of eigenvalue j
    b 		= imag(lambda(j)); 	% imaginary part of eigenvalue j
    abs_lambda_sq= abs(lambda(j))^2;  	% squared absolute value of eigenvalue j
    tau(1,j) 	= -2/log(abs_lambda_sq);% damping time of eigenmode j

    % Excitation of eigenmode j 
    exctn(j) 	= real(CovDcpld(j,j) / (1-abs_lambda_sq)); 

    % Assemble derivative of eigenvalue with respect to parameters in 
    % the coefficient matrix A 
    dot_lam 	= zeros(m^2*p, 1);
    for k=1:m
      dot_lam(k:m:k+(m*p-1)*m) = BigS_inv(j,k) .* BigS(1:m*p,j);
    end
    dot_a 	= real(dot_lam); 	% derivative of real part of lambda(j)
    dot_b 	= imag(dot_lam); 	% derivative of imag part of lambda(j)
    
    % Derivative of the damping time tau w.r.t. parameters in A
    phi 	= tau(1,j)^2 / abs_lambda_sq * (a*dot_a + b*dot_b);
    % Margin of error for damping time tau
    tau(2,j) 	= t * sqrt(phi'*Sigma_A*phi);
        
    % Period of eigenmode j and margin of error for period. (The
    % if-statements avoid warning messages that may otherwise result
    % from a division by zero)
    if (b == 0 & a >= 0)      % purely real, nonnegative eigenvalue
      per(1,j)	= Inf;   		
      per(2,j)  = 0;         
    elseif (b == 0 & a < 0)   % purely real, negative eigenvalue
      per(1,j)	= 2;     		
      per(2,j)  = 0;         
    else                      % complex eigenvalue
      per(1,j)	= 2*pi/abs(atan2(b,a)); 
      
      % Derivative of period with respect to parameters in A
      phi 	= per(1,j)^2 / (2*pi*abs_lambda_sq)*(b*dot_a-a*dot_b);
      % Margin of error for period
      per(2,j) 	= t * sqrt(phi'*Sigma_A*phi);
    end
  end

  % Give the excitation as `relative importance' that sums to one
  exctn 	= exctn/sum(exctn);
  
  % Compute confidence intervals for eigenmodes 
  % -------------------------------------------
  % Shorthands for matrix products
  XX 	        = real(BigS)'*real(BigS);
  YY 	        = imag(BigS)'*imag(BigS);
  XY 	        = real(BigS)'*imag(BigS);

  % Need confidence intervals only for last m rows of BigS
  row1 	        = (p-1)*m+1; % first row for which confidence interval is needed

  mp = m*p;                           	% dimension of equivalent AR(1) model
  for k=1:mp        			% loop over columns of S
    for l=row1:mp  			% loop over rows of S
					
      % evaluate gradient of S_{lk}
      for ii=1:m
	for jj=1:mp
	  % compute derivative with respect to A(ii,jj)
	  zsum = 0;
	  zkkr = 0; 			% real part of Z_{kk}
	  zkki = 0; 			% imaginary part of Z_{kk}
	  for j=1:mp
	    if (j ~= k) 		% sum up elements that appear in Z_{kk}
	      zjk  = BigS_inv(j,ii)*BigS(jj,k)/(lambda(k)-lambda(j));
	      zjkr = real(zjk);
	      zjki = imag(zjk);
	      zkkr = zkkr+zjki*(XY(k,j)-XY(j,k))-zjkr*(XX(k,j)+YY(k,j));
	      zkki = zkki+zjki*(YY(k,j)-XX(k,j))-zjkr*(XY(k,j)+XY(j,k));
	      zsum = zsum+BigS(l,j)*zjk;
	    end
	  end
	  % now add Z_{kk}
	  zkki = zkki / (XX(k,k)-YY(k,k));
	  grad_S((jj-1)*m+ii) = zsum+BigS(l,k)*(zkkr+i*zkki);
	end     
      end 
      Serr(l,k) = t * ( sqrt( real(grad_S)*Sigma_A*real(grad_S)') ...
			 + i*sqrt( imag(grad_S)*Sigma_A*imag(grad_S)') );
    end
  end

  % Only return last m*p rows of Serr
  Serr = Serr(row1:m*p, :);






