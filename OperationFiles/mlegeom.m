% n=1 : mean of mle gaussian fit

function out=mlegeom(y,n)
phat=mle(y,'distribution','geometric');
out=phar(n);