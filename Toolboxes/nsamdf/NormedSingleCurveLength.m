function [rx1,pxy,frq] = NormedSingleCurveLength(x,lag,fs,figurep,nrmdegree)

% this function computes crossed curve length for single time series x

rx1 = zeros(1,lag+1);

for delay=0:lag
    rx1(delay+1) = norm(x(1:end-delay) - x(delay+1:end),nrmdegree);
end

rx = [fliplr(rx1) rx1(2:end)];

pxy = abs(fft(rx));

frq = linspace(0,0.5*fs,lag+1);

if figurep == 1
    %figure, plot(rx,'k')
    [~,ind] = min( abs(frq - 60) );

    figure, plot(frq(2:ind), pxy(2:ind),'b')
end

end
