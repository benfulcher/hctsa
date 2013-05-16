function out=CO_fzcac(y)

tau = 0;
% AC1=CO_autocorr(y,1);
AC1 = 1;
while tau<min(400,length(y)-1);
    tau = tau+1;
    AC2 = CO_autocorr(y,tau);
    if AC1*AC2<0 % zero crossing
        out = tau;
        return
    else
        AC1 = AC2;
    end
end

out = tau;

end