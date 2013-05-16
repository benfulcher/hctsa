function out = CO_ami_noise(y,tau,meth,nbins)
% Looks at how much AMI curves change as noise is added
% Ben Fulcher September 2009
% y needs to be z-scored

noiser = (0:0.1:2);
nr = length(noiser);
amis = zeros(nr,1);
noise = randn(size(y));

for i = 1:nr;
    amis(i) = CO_ami_benhist(y+noiser(i)*noise,tau,meth,nbins);
end

% plot(noiser,amis,'.-k');
% out=amis;

out.pdec = length(find(diff(amis)<0))/(nr-1);
out.meanch = mean(diff(amis));
out.ac1 = CO_autocorr(amis,1);
out.ac2 = CO_autocorr(amis,2);

% fit exponential decay to output
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[amis(1) -1]);
f = fittype('a*exp(b*x)','options',s);
[c,gof] = fit(noiser',amis,f);

% plot
% cc = bengetcmap('set1',2,1);
% % figure('color','w');
% hold on; box('on')
% plot(noiser,c.a*exp(c.b*noiser),'color',cc{2},'linewidth',2)
% plot(noiser,amis,'.-','color',cc{1})
% xlabel('\eta'); ylabel('AMI_1')

out.fitexpa = c.a;
out.fitexpb = c.b;
out.fitexpr2 = gof.rsquare;
out.fitexpadjr2 = gof.adjrsquare;
out.fitexprmse = gof.rmse;

% fit linear to output
p = polyfit(noiser',amis,1);
out.fitlina = p(1);
out.fitlinb = p(2);
linfit = polyval(p,noiser);
out.mse = mean((linfit' - amis).^2);

% crosses mean
out.pcrossmean = length(sgnchange(amis-mean(amis)))/(nr-1);

end