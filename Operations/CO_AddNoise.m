% CO_AddNoise
% 
% Analyzes changes in the automutual information function with the addition of
% noise to the input time series.
% Adds Gaussian-distributed noise to the time series with increasing standard
% deviation, eta, across the range eta = 0, 0.1, ..., 2, and measures the
% mutual information at each point using histograms with
% nbins bins (implemented using CO_ami_benhist).
% 
% The output is a set of statistics on the resulting set of automutual
% information estimates, including a fit to an exponential decay, since the
% automutual information decreases with the added white noise.
% 
% Can calculate these statistics for time delays 'tau', and for a number 'nbins'
% bins.
% 
% This algorithm is quite different, but was based on the idea of 'noise
% titration' presented in: "Titration of chaos with added noise", Chi-Sang Poon
% and Mauricio Barahona P. Natl. Acad. Sci. USA, 98(13) 7107 (2001)
% 

function out = CO_AddNoise(y,tau,meth,nbins)
% Ben Fulcher, September 2009

%% Check inputs
if nargin < 2
    tau = []; % set default in CO_ami_benhist
end
if nargin < 3
    meth = ''; % set default in CO_ami_benhist
end
if nargin < 4
    nbins = [];
end

% Preliminaries
noiser = (0:0.1:2); % across this noise range
nr = length(noiser);
amis = zeros(nr,1);
noise = randn(size(y));

for i = 1:nr
    amis(i) = CO_ami_benhist(y+noiser(i)*noise,tau,meth,nbins);
end

% plot(noiser,amis,'.-k');
% out = amis;

out.pdec = sum(diff(amis) < 0)/(nr-1);
out.meanch = mean(diff(amis));
out.ac1 = CO_autocorr(amis,1);
out.ac2 = CO_autocorr(amis,2);

% Fit exponential decay to output
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[amis(1) -1]);
f = fittype('a*exp(b*x)','options',s);
[c, gof] = fit(noiser',amis,f);

% plot
% cc = BF_getcmap('set1',2,1);
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

% Fit linear function to output
p = polyfit(noiser',amis,1);
out.fitlina = p(1);
out.fitlinb = p(2);
linfit = polyval(p,noiser);
out.mse = mean((linfit' - amis).^2);

% Number of times the AMI function crosses its mean
out.pcrossmean = sum(BF_sgnchange(amis-mean(amis)))/(nr-1);

end