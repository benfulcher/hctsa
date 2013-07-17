function out = SY_dynstdder(y,maxd)
% Operation inspired by a comment on the Matlab Central forum:
% "You can measure the standard deviation of the n-th derivative, if you
% like." -- Vladimir Vassilevsky, DSP and Mixed Signal Design Consultant
% from http://www.mathworks.de/matlabcentral/newsreader/view_thread/136539

doplot = 0;

if nargin < 2 || isempty(maxd)
    maxd = 10; % do 10 by default
end

ms = zeros(maxd,1);
for i = 1:maxd
    ms(i) = SY_stdnthder(y,i);
end

if doplot
    plot(ms)
end

% fit exponential growth/decay
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1, 0.5*sign(ms(end)-ms(1))]);
f = fittype('a*exp(b*x)','options',s);
[c, gof] = fit((1:maxd)',ms,f);
out.fexp_a = c.a;
out.fexp_b = c.b; % this is important
out.fexp_r2 = gof.rsquare; % this is more important!
out.fexp_adjr2 = gof.adjrsquare;
out.fexp_rmse = gof.rmse;

% looks like excellent fit to exponential: regular signals decrease, irregular signals
% increase...

end