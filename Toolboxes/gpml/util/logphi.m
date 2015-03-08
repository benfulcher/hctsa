% Safe computation of logphi(z) = log(normcdf(z)) and its derivatives
%                    dlogphi(z) = normpdf(x)/normcdf(x).
% The function is based on index 5725 in Hart et al. and gsl_sf_log_erfc_e.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2013-11-13.
function [lp,dlp,d2lp,d3lp] = logphi(z)
  z = real(z);                                 % support for real arguments only
  lp = zeros(size(z));                                         % allocate memory
  id1 = z.*z<0.0492;                                 % first case: close to zero
  lp0 = -z(id1)/sqrt(2*pi);
  c = [ 0.00048204; -0.00142906; 0.0013200243174; 0.0009461589032;
       -0.0045563339802; 0.00556964649138; 0.00125993961762116;
       -0.01621575378835404; 0.02629651521057465; -0.001829764677455021;
       2*(1-pi/3); (4-pi)/3; 1; 1];
  f = 0; for i=1:14, f = lp0.*(c(i)+f); end, lp(id1) = -2*f-log(2);
  id2 = z<-11.3137;                                    % second case: very small
  r = [ 1.2753666447299659525; 5.019049726784267463450;
        6.1602098531096305441; 7.409740605964741794425;
        2.9788656263939928886 ];
  q = [ 2.260528520767326969592;  9.3960340162350541504;
       12.048951927855129036034; 17.081440747466004316; 
        9.608965327192787870698;  3.3690752069827527677 ];
  num = 0.5641895835477550741; for i=1:5, num = -z(id2).*num/sqrt(2) + r(i); end
  den = 1.0;                   for i=1:6, den = -z(id2).*den/sqrt(2) + q(i); end
  e = num./den; lp(id2) = log(e/2) - z(id2).^2/2;
  id3 = ~id2 & ~id1; lp(id3) = log(erfc(-z(id3)/sqrt(2))/2);  % third case: rest
  if nargout>1                                        % compute first derivative
    dlp = zeros(size(z));                                      % allocate memory
    dlp( id2) = abs(den./num) * sqrt(2/pi); % strictly positive first derivative
    dlp(~id2) = exp(-z(~id2).*z(~id2)/2-lp(~id2))/sqrt(2*pi); % safe computation
    if nargout>2                                     % compute second derivative
      d2lp = -dlp.*abs(z+dlp);             % strictly negative second derivative
      if nargout>3                                    % compute third derivative
        d3lp = -d2lp.*abs(z+2*dlp)-dlp;     % strictly positive third derivative
      end
    end
  end