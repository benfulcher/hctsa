
function out=CO_fmac(y)

a=[];
for i=0:length(y)-1
    a(i+1)=CO_autocorr(y,i);
    if i>1 && a(i-1)-a(i)>0 && a(i)-a(i+1)<0; % minimum
        out=i-1;
        return
    end
end
% if no minimum in entire spectrum
out=length(y);
end