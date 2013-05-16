function ox=adjph(x)
%ADJPH	Normalization of columns of a complex matrix.
%
%  Given a complex matrix X, OX=ADJPH(X) returns the complex matrix OX
%  that is obtained from X by multiplying column vectors of X with
%  phase factors exp(i*phi) such that the real part and the imaginary
%  part of each column vector of OX are orthogonal and the norm of the
%  real part is greater than or equal to the norm of the imaginary
%  part.
%
%  ADJPH is called by ARMODE.
%
%  See also ARMODE.

%  Modified 16-Dec-99
%  Author: Tapio Schneider
%	   tapio@gps.caltech.edu

  for j = 1:size(x,2)				
    a       = real(x(:,j));                     % real part of jth column of x
    b       = imag(x(:,j));                     % imag part of jth column of x
    phi     = .5*atan( 2*sum(a.*b)/(b'*b-a'*a) );
    bnorm   = norm(sin(phi).*a+cos(phi).*b);    % norm of new imaginary part
    anorm   = norm(cos(phi).*a-sin(phi).*b);    % norm of new real part
    if bnorm > anorm 
      if phi < 0
	phi = phi-pi/2;
      else
	phi = phi+pi/2;
      end
    end
    ox(:,j) = x(:,j).*exp(i*phi);
  end






