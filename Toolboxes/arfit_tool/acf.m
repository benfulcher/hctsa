function []=acf(x, k, caption)
%ACF	Plot of sample autocorrelation function.
%
%  ACF(x) plots the sample autocorrelation function of the univariate
%  time series in vector x.  By default, the sample autocorrelation
%  function is plotted up to lag 25. ACF(x,k) plots the sample
%  autocorrelation function up to lag k. ACF(X,k,'name') sets the
%  title of the plot to 'name'.
%
%  The approximate 95% confidence limits of the autocorrelation
%  function of an IID process of the same length as X are also
%  displayed.  Sample autocorrlations lying outside the 95% confidence
%  intervals of an IID process are marked by an asterisk.
%
%  ACF requires XCORR from the Signal Processing Toolbox.
%
%  See also XCORR.

%  Modified 30-Dec-99
%  Author: Tapio Schneider
%	   tapio@gps.caltech.edu

  if ~exist('xcorr')
    error('ACF requires XCORR from the Signal Processing Toolbox.')
  end   

  [m,n]	  = size(x);
  if (min(m,n) > 1) error('Time series must be univariate.'); end
  n 	  = max(m,n);

  if (nargin < 3) caption='ACF'; end
  if (nargin < 2) k=25; end
  if (nargin < 1) error('Usage: ACF(vector).'); end

  % Compute autocorrelation sequence with XCORR from the Signal 
  % Processing Toolbox
  cor     = xcorr(x,'coeff');
  rho     = cor(n:n+k);       % autocorrelation function up to lag k
  abscis  = [0:1:k];          % abscissa for plot

  % Approximate 95% confidence limits for IID process of the same length
  bound   = zeros(size(rho));
  bound   = bound + 1.96/sqrt(n);

  % Initialize abscissas and ordinates for ACF values within and outside the
  % approximate 95% confidence intervals of IID process
  inabsc  = zeros(1, k+1);
  inl     = zeros(1, k+1);
  outabsc = zeros(1, k+1);
  outl    = zeros(1, k+1);

  % Find lags within and outside approximate 95% confidence
  % intervals; start with lag 0
  inl(1)  = rho(1);           % stores ACF at lags within cfd intervals
  inabsc(1) = 0;              % lags at which ACF is within cfd intervals
  nin     = 1;                % number of points within confidence intervals 
  nout    = 0;                % number of points outside confidence intervals 
  for i=2:k+1 
    if abs(rho(i)) > bound(1) % point outside confidence intervals
      nout          = nout+1;
      outl(nout)    = rho(i);
      outabsc(nout) = i-1;
    else                      % point within confidence intervals
      nin           = nin+1;
      inl(nin)      = rho(i);
      inabsc(nin)   = i-1;
    end;
  end;

  % Plot ACF
  plot(abscis, rho, abscis, bound, '-.', abscis, -bound, '-.', ...
       outabsc(1:nout), outl(1:nout), '*', ...
       inabsc(1:nin),   inl(1:nin),   'o');
  axis([0 k -1 1])
  title(caption)
  xlabel('Lag')
  ylabel('Autocorrelation function')









