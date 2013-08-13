% function e = MS_rms(y);
%
% e is the l2-norm of row vector y, for a n-by-m matrix e is the n-by-1 column 
% vector which is the l2-norm of the n rows of y.;
%
% Michael Small
% michael.small@uwa.edu.au, http://school.maths.uwa.edu.au/~small/
% 3/3/2005
% For further details, please see M. Small. Applied Nonlinear Time Series
% Analysis: Applications in Physics, Physiology and Finance. Nonlinear Science
% Series A, vol. 52. World Scientific, 2005. (ISBN 981-256-117-X) and the
% references therein.
% (Minor edits by Ben Fulcher, 2010)

function e = MS_rms(y);

if ndims(y) > 2
  
  e = sqrt(mean(y.^2));
  
else

  [a, b] = size(y);

  if b == 1
    e = abs(y);
  else
    if a == 1
      y = y';
    end;
    e = sqrt(mean(y'.^2))';
  end
end
end