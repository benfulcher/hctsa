
function out = CO_f1ecac(y)
oone = 1/exp(1);
a(1) = CO_information(y,y); % very weird -- why is this not 1?? Or use autocor?

for i = 2:length(y)-1
    a(i) = CO_autocorr(y,i);
    if (a(i-1)-oone)*(a(i)-oone)<0
        out=i;
        return
    end
end

% if no minimum in entire spectrum
out = length(y);

end

% 
% tau=1;
% 
% AC1=CO_autocorr(y,1);
% while tau<min(400,length(y)-1);
%     tau=tau+1;
%     AC2=CO_autocorr(y,tau);
%     if (AC1-oone)*(AC2-oone)<0 % zero crossing
%         out=tau;
%         return
%     else
%         AC1=AC2;
%     end
% end
% out=tau;
% end