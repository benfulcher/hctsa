function out = EN_dynhisten(y,ksorhist)
% Ranges over bin size and returns statistics on the output
% Ben Fulcher, 2009

doplot = 0; % plot outputs

if nargin < 2
    ksorhist = 'hist'; % use histogram by default
end

switch ksorhist
    case 'ks'
        widthr = (0.01:0.01:1);
        hs = zeros(length(widthr),1);
        for i = 1:length(widthr)
            width = widthr(i);
            [px, xr] = ksdensity(y,'width',width,'function','pdf');
            hs(i) = -sum(px.*log(eps+px))*(xr(2)-xr(1));
        end
        if doplot
            plot(widthr,hs);
        end
    case 'hist'
        binsizer = (2:50); % range of binsizes
        hs = zeros(length(binsizer),1);
        for i = 1:length(binsizer)
            nbins = binsizer(i);
            [px, xr] = hist(y,nbins);
            px = px/(sum(px)*(xr(2)-xr(1)));  
            hs(i) = -sum(px.*log(eps+px))*(xr(2)-xr(1));
        end
        if doplot
            plot(binsizer,hs);
        end
otherwise
    error('Unknown distribution method: %s',ksorhist)
end

% Summary statistics
out.ch = hs(end)-hs(1);
out.pch = (hs(end)-hs(1))/abs(hs(1));
out.mean = mean(hs);
out.std = std(hs);

% Fit Exponential decay to absolute ACF using nonlinear least squares:
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[2, -0.5]);
f = fittype('a*exp(b*x)','options',s);
[c, gof] = fit(binsizer',hs,f);

out.fexpabsacf_a = c.a;
out.fexpabsacf_b = c.b; % this is important
out.fexpabsacf_r2 = gof.rsquare; % this is more important!
out.fexpabsacf_adjr2 = gof.adjrsquare;
out.fexpabsacf_rmse = gof.rmse;

% (**) look at variance of residuals

end