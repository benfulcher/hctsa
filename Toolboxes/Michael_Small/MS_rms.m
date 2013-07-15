function e = MS_rms(y);

% function e=MS_rms(y);
%
% e is the l2-norm of row vector y, for a n-by-m matrix e is the n-by-1 column 
% vector which is the l2-norm of the n rows of y.;
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk

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