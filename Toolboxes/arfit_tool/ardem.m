%ARDEM	Demonstrates modules of the ARfit package.

%   Revised: 30-Dec-99 Tapio Schneider 

format short
format compact

echo on
clc

%  ARfit is a collection of Matlab modules for the modeling of
%  multivariate time series with autoregressive (AR) models. ARfit
%  contains modules for estimating parameters of AR models from given
%  time series data; for checking the adequacy of an estimated
%  AR model; for analyzing eigenmodes of an estimated AR model; and
%  for simulating AR processes.
%
%  This demo illustrates the use of ARfit with a bivariate AR(2)
%  process
%
%       v(k,:)' = w' + A1*v(k-1,:)' + A2*v(k-2,:)' + eta(k,:)',
%
%  where the row vectors eta(k,:) are independent and identically
%  distributed Gaussian noise vectors with zero mean and covariance
%  matrix C. The kth row v(k,:) of the 2-column matrix v represents an
%  observation of the process at instant k. The intercept vector w is
%  included to allow for a nonzero mean of the AR(p) process.

pause   	% Press any key to continue.

clc

%  Let us simulate observations from a bivariate AR(2) process,
%  choosing the parameters

w = [ 0.25 ;  0.1 ];

% for the intercept vector,

A1 = [ 0.4   1.2;   0.3   0.7 ];

%  and

A2 = [ 0.35 -0.3;  -0.4  -0.5 ];

%  for the AR coefficient matrices, and 

C = [ 1.00  0.50;   0.50  1.50 ];

%  for the noise covariance matrix. The two 2x2 matrices A1 and A2 are
%  assembled into a single 2x4 coefficient matrix:

A = [ A1 A2 ];

%  We use the module ARSIM to simulate 200 observations of this AR
%  process:

v = arsim(w, A, C, 200);

pause   	% Press any key to continue.

clc

%  Suppose that we have no information about how the time series
%  data v are generated, but we want to try to fit an AR model to the
%  time series. That is, we must estimate the AR model parameters,
%  including the model order. Assuming that the correct model order
%  lies between 1 and 5, we use the module ARFIT to determine the
%  optimum model order using Schwarz's Bayesian Criterion (SBC):

pmin = 1;
pmax = 5;
[west, Aest, Cest, SBC, FPE, th] = arfit(v, pmin, pmax);

%  The output arguments west, Aest, and Cest of ARFIT are the
%  estimates of the intercept vector w, of the coefficient matrix A,
%  and of the noise covariance matrix C. (The matrix th will be needed
%  later in the computation of confidence intervals.) These parameters
%  are estimated for a model of order popt, where the optimum model
%  order popt is chosen as the minimizer of an approximation to
%  Schwarz's Bayesian Criterion. The selected order popt in our
%  example is:

m    = 2;                 % dimension of the state space
popt = size(Aest, 2) / m;

disp(['popt = ', num2str(popt)])

pause   	% Press any key to continue.

clc

%  Besides the parameter estimates for the selected model, ARFIT has
%  also returned approximations to Schwarz's Bayesian Criterion SBC
%  and to Akaike's Final Prediction Error FPE, each for models of
%  order pmin:pmax. In this demo, the model order was chosen as the
%  minimizer of SBC. Here are the SBCs for the fitted models of order
%  1,...,5:

disp(SBC)
 
%  To see if using Akaike's Final Prediction Error as a criterion to
%  select the model order would have resulted in the same optimum
%  model order, compare the FPEs:

disp(FPE) 

%  Employing FPE as the order selection criterion, the optimum model
%  order would have been chosen as the minimizer of FPE. The values of
%  the order selection criteria are approximations in that in
%  evaluating an order selection criterion for a model of order p <
%  pmax, the first pmax-p observations of the time series are ignored.

pause   	% Press any key to continue.

clc

%  Next it is necessary to check whether the fitted model is adequate
%  to represent the time series v. A necessary condition for model
%  adequacy is the uncorrelatedness of the residuals. The module ARRES
 
[siglev,res] = arres(west,Aest,v);

%  returns the time series of residuals res as well as the
%  significance level siglev with which a modified Li-McLeod
%  portmanteau test rejects the null hypothesis that the residuals are
%  uncorrelated. A model passes this test if, say, siglev > 0.05.  In
%  our example, the significance level of the modified Li-McLeod
%  portmanteau statistic is

disp(siglev); 

%  (If siglev is persistently smaller than about 0.05, even over
%  several runs of this demo, then there is most likely a problem with
%  the random number generator that ARSIM used in the simulation of
%  v.)

pause   	% Press any key to continue.

clc

if exist('xcorr')    % If the Signal Processing Toolbox is installed, plot
                     % autocorrelation function of residuals ...

%  Using ACF, one can also plot the autocorrelation function of, say, 
%  the first component of the time series of residuals:

acf(res(:,1));

%  95% of the autocorrelations for lag > 0 should lie between the
%  dashdotted confidence limits for the autocorrelations of an IID
%  process of the same length as the time series of residuals res.

%  Since uncorrelatedness of the residuals is only a necessary
%  condition for model adequacy, further diagnostic tests should, in
%  practice, be performed. However, we shall go on to the estimation
%  of confidence intervals.

pause   	% Press any key to continue.
end

clc

%  Being reasonably confident that the model adequately represents the
%  data, we compute confidence intervals for the AR parameters with
%  the module ARCONF:

[Aerr, werr] = arconf(Aest, Cest, west, th);

%  The output arguments Aerr and werr contain margins of error such
%  that (Aest +/- Aerr) and (west +/- werr) are approximate 95%
%  confidence intervals for the individual AR coefficients and for the
%  components of the intercept vector w. Here is the estimated
%  intercept vector with margins of error for the individual
%  components in the second column:

disp([west werr])

%  For comparison, the `true' vector as used in the simulation:

disp(w)

pause   	% Press any key to display the other parameter estimates.

clc

%  The estimated coefficient matrix:

disp(Aest)

%  with margins of error for its elements:

disp(Aerr)

%  For comparison, the `true' AR coefficients:

disp(A)

pause   	% Press any key to continue.
		
echo off
%  Compute `true' eigenmodes from model parameters:
%  Eigenvectors and eigenvalues of corresponding AR(1) coefficient
%  matrix: 
[Strue,EvTrue] = eig([A; eye(2) zeros(2)]);  
EvTrue         = diag(EvTrue)';   % `true' eigenvalues
Strue          = adjph(Strue);    % adjust phase of eigenmodes

clc
echo on

%  Finally, ARMODE computes the eigendecomposition of the fitted AR
%  model:

[S, Serr, per, tau, exctn] = armode(Aest, Cest, th); 

%  The columns of S are the estimated eigenmodes:

disp(S)

%  with margins of error Serr:

disp(Serr)

%  The intervals (S +/- Serr) are approximate 95% confidence intervals
%  for the individual components of the eigenmodes.  Compare the
%  estimated eigenmodes above with the eigenmodes obtained from the
%  `true' model parameters:

disp(Strue(3:4,:))

%  (Note that the estimated modes can be a permutation of the `true'
%  modes. The sign of the modes is also ambiguous.)

pause   	% Press any key to continue.

echo off
pertrue = 2*pi./abs(angle(EvTrue)); % `true' periods

clc

echo on
%  Associated with the eigenmodes are the following oscillation periods:

disp(per) 

%  The second row contains margins of error for the periods such that
%  (per(1,k) +/- per(2,k)) is an approximate 95% confidence interval
%  for the period of eigenmode S(:,k). [Note that for a purely
%  relaxatory eigenmode, the period is infinite (Inf).] Compare the
%  above estimated periods with the `true' periods:

disp(pertrue)

pause   	% Press any key to get the damping time scales.

echo off
tautrue = -1./log(abs(EvTrue)); % `true' damping time scales

clc

echo on 

%  The damping times associated with each eigenmode are:

disp(tau)

%  with margins of error again in the second row. For comparison, the
%  `true' damping times:

disp(tautrue)

pause   	% Press any key to get the excitation of each eigenmode.
echo off
%  Compute `true' excitation of eigenmodes from the designed parameters:
p  = 2;              % true model order

invStr = inv(Strue); % inverse of matrix with eigenvectors as columns

% covariance matrix of corresponding decoupled AR(1) system
CovDcpld = invStr*[C zeros(2,(p-1)*2); zeros((p-1)*2, p*2)]*invStr';

% diagonal of that covariance matrix
DgCovDcpld = real(diag(CovDcpld))';

% excitation 
TrueExctn = DgCovDcpld(1:2*p)./(1-abs(EvTrue).^2);
 
% normalize excitation 
TrueExctn = TrueExctn./sum(TrueExctn);            

clc

echo on
%  ARMODE has also returned the excitations, measures of the relative
%  dynamical importance of the eigenmodes:

disp(exctn)

%  Compare the estimated excitations with the `true' excitations
%  computed from the parameters used in the simulation:

disp(TrueExctn)

echo off
disp('End')

