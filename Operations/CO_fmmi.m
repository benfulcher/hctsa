function out = CO_fmmi(y)
% CO_MI
% Ben Fulcher in early days 2008 -- outsources CO_information code that 
% may be shit...

a = [];
for i = 0:length(y)-1
    a(i+1) = CO_information(y(1:end-i),y(1+i:end));
    if i>1 && a(i-1)-a(i)>0 && a(i)-a(i+1)<0; % minimum
        out = i;
        return
    end
end
% if no minimum in entire spectrum
out = length(y);

end