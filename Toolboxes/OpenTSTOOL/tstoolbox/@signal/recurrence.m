function rs=recurrence(s,func)

% tstoolbox/signal/recurrence
%

d=data(s);
[N,D]=size(d);


result=zeros(N);
for i=1:N
  for i1=1:N
    a=d(i,:);
    b=d(i1,:);
    result(i,i1)=eval(func);
  end
end
rs=signal(result);
rs=setplothint(rs,'spectrogram');

return
