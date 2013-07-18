function tau=MS_firstzero(y);
  
%function tau=firstzero(y);
%
%find the first zero of the autocorrelation function of y.
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk
  
len=50;  
[r,t]=acorr(y,len);
lY=length(y);
%lY=1000;

while (all(r(2:(end-1))>0)),
  len=len*2;
  if (len==2*lY),
    disp('WARNING : No minium found in firstzero(y)');
    tau=0;
    return;
  end;
  if (len>lY)
    len=lY;
  end;
  [r,t]=acorr(y,len);
end;

r=r((len+1):end);

t0=find(r>0);
t0=t0(find(diff(t0)>1));

if length(t0)>0,
  t0=t0(1);
  t1=t0+1;
else
  t0=min(find(r<0));
  t1=t0+1;
end;

if t0==(len+1),
  tau=t0;
  return;
end;

if (-r(t1)>=r(t0))
  tau=t0;
else,
  tau=t1;
end;
tau=tau-1;

%if tau>(len/4),
%  tau=1;
%end;
