function out=sgnchange(y)

out=find(y(2:end).*y(1:end-1)<0);

end