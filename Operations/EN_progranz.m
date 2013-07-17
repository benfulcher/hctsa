function out = EN_progranz(y,howtorand)
% Progressively randomizes the input time series, returning measures of
% what happens as a result
% Inputs:
%        y: the (zscored) time series
%        howtorand:
%             1) 'statdist' -- substitutes a random element of the time
%        series with one from the original time series distribution
%             2) 'dyndist' -- overwrites a random element of the time
%             series with another random element
%             3) 'permute' -- progressively permutes elements of the time
%             series
% Ben Fulcher, October 2009

if ~BF_iszscored(y)
    warning('The input time series should be z-scored for EN_progranz')
end

if nargin < 2 || isempty(howtorand)
    howtorand = 'statdist'; % use statdist by default
end

N = length(y); % length of the time series

nstats = 11;
randp_max = 2; % time series has been randomized to 2 times its length
rand_inc = 10/100; % this proportion of the time series between calculations
ncalcs = randp_max/rand_inc;
calc_ints = floor(randp_max*N/ncalcs);
calc_pts = (0:calc_ints:randp_max*N);
ncalcs = length(calc_pts); % some rounding issues inevitable

stats = zeros(ncalcs,nstats); % record a stat at each randomization increment

y_rand = y; % this vector will be randomized

stats(1,:) = doyourcalcthing(y,y_rand); % initial condition: apply on itself

for i = 1:N*randp_max
    switch howtorand
        case 'statdist'
            % randomize by substituting a random element of the time series by
            % a random element from the static original time series distribution
            y_rand(randi(N)) = y(randi(N));
        case 'dyndist'
            % randomize by substituting a random element of the time series
            % by a random element of the current, already partially randomized, 
            % time series
            y_rand(randi(N)) = y_rand(randi(N));
        case 'permute'
            % randomize by permuting elements of the time series so that
            % the distribution remains static
            randis = randi(N,[2, 1]);
            y_rand(randis(1)) = y_rand(randis(2));
            y_rand(randis(2)) = y_rand(randis(1));
    end
    
    if ismember(i,calc_pts)
        stats(calc_pts == i,:) = doyourcalcthing(y,y_rand);
        % disp([num2str(i) ' out of ' num2str(N*randp_max)])
    end
    
end

% plot(stats,'.-');
% out=stats;

%% Fit to distributions of output
r = (1:size(stats,1))'; % gives an 'x-axis' for fitting

% 1) xcn1: cross correlation at negative 1
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[stats(1,1) -0.1]);
f = fittype('a*exp(b*x)','options',s);
[c, gof] = fit(r,stats(:,1),f);
out.xcn1fexpa = c.a;
out.xcn1fexpb = c.b;
% out.xcn1fexpc = c.c;
out.xcn1fexpr2 = gof.rsquare;
out.xcn1fexpadjr2 = gof.adjrsquare;
out.xcn1fexprmse = gof.rmse;

out.xcn1diff = abs((stats(end,1)-stats(1,1))/stats(end,1));
out.xcn1hp = gethp(stats(:,1));

% 2) xc1: cross correlation at lag 1
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[stats(1,2) -0.1]);
f = fittype('a*exp(b*x)','options',s);
[c, gof] = fit(r,stats(:,2),f);
out.xc1fexpa = c.a;
out.xc1fexpb = c.b;
% out.xc1fexpc = c.c;
out.xc1fexpr2 = gof.rsquare;
out.xc1fexpadjr2 = gof.adjrsquare;
out.xc1fexprmse = gof.rmse;

out.xc1diff = abs((stats(end,2)-stats(1,2))/stats(end,2));
out.xc1hp = gethp(stats(:,2));

% 3) d1: norm of differences between original and randomized time series
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[-stats(end,3) -0.2 stats(end,3)]);
f = fittype('a*exp(b*x)+c','options',s);
[c, gof] = fit(r,stats(:,3),f);
out.d1fexpa = c.a;
out.d1fexpb = c.b;
out.d1fexpc = c.c;
out.d1fexpr2 = gof.rsquare;
out.d1fexpadjr2 = gof.adjrsquare;
out.d1fexprmse = gof.rmse;

out.d1diff = abs((stats(end,3)-stats(1,3))/stats(end,3));
out.d1hp = gethp(stats(:,3));

% 4) ac1
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[stats(1,4) -0.2]);
f = fittype('a*exp(b*x)','options',s);
[c, gof] = fit(r,stats(:,4),f);
out.ac1fexpa = c.a;
out.ac1fexpb = c.b;
% out.ac1fexpc = c.c;
out.ac1fexpr2 = gof.rsquare;
out.ac1fexpadjr2 = gof.adjrsquare;
out.ac1fexprmse = gof.rmse;

out.ac1diff = abs((stats(end,4)-stats(1,4))/stats(end,4));
out.ac1hp = gethp(stats(:,4));


% 5) ac2
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[stats(1,5) -0.2]);
f = fittype('a*exp(b*x)','options',s);
[c, gof] = fit(r,stats(:,5),f);
out.ac2fexpa = c.a;
out.ac2fexpb = c.b;
% out.ac2fexpc = c.c;
out.ac2fexpr2 = gof.rsquare;
out.ac2fexpadjr2 = gof.adjrsquare;
out.ac2fexprmse = gof.rmse;

out.ac2diff = abs((stats(end,5)-stats(1,5))/stats(end,5));
out.ac2hp = gethp(stats(:,5));

% 6) ac3
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[stats(1,6) -0.2]);
f = fittype('a*exp(b*x)','options',s);
[c, gof] = fit(r,stats(:,6),f);
out.ac3fexpa = c.a;
out.ac3fexpb = c.b;
% out.ac3fexpc = c.c;
out.ac3fexpr2 = gof.rsquare;
out.ac3fexpadjr2 = gof.adjrsquare;
out.ac3fexprmse = gof.rmse;

out.ac3diff = abs((stats(end,6)-stats(1,6))/stats(end,6));
out.ac3hp = gethp(stats(:,6));

% 7) ac4
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[stats(1,7) -0.4]);
f = fittype('a*exp(b*x)','options',s);
[c, gof] = fit(r,stats(:,7),f);
out.ac4fexpa = c.a;
out.ac4fexpb = c.b;
% out.ac4fexpc = c.c;
out.ac4fexpr2 = gof.rsquare;
out.ac4fexpadjr2 = gof.adjrsquare;
out.ac4fexprmse = gof.rmse;

out.ac4diff = abs((stats(end,7)-stats(1,7))/stats(end,7));
out.ac4hp = gethp(stats(:,7));

% 8) shen (Shannon entropy)
% I think this is all rubbish
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1, -1, 0.4]);
f = fittype('a*exp(b*x)+c','options',s);
[c, gof] = fit(r,stats(:,8),f);
out.shenfexpa = c.a;
out.shenfexpb = c.b;
out.shenfexpc = c.c;
out.shenfexpr2 = gof.rsquare;
out.shenfexpadjr2 = gof.adjrsquare;
out.shenfexprmse = gof.rmse;

out.shendiff = abs((stats(end,8)-stats(1,8))/stats(end,8));
out.shenhp = gethp(stats(:,8));

% 9) sampen (Sample Entropy)
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[-stats(end,9) -0.2 stats(end,9)]);
f = fittype('a*exp(b*x)+c','options',s);
[c, gof] = fit(r,stats(:,9),f);
out.sampen2_02fexpa = c.a;
out.sampen2_02fexpb = c.b;
out.sampen2_02fexpc = c.c;
out.sampen2_02fexpr2 = gof.rsquare;
out.sampen2_02fexpadjr2 = gof.adjrsquare;
out.sampen2_02fexprmse = gof.rmse;

out.sampen2_02diff = abs((stats(end,9)-stats(1,9))/stats(end,9));
out.sampen2_02hp = gethp(stats(:,9));

% 10) statav5
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[-stats(end,10) -0.1 stats(end,10)]);
f = fittype('a*exp(b*x)+c','options',s);
[c, gof] = fit(r,stats(:,10),f);
out.statav5fexpa = c.a;
out.statav5fexpb = c.b;
out.statav5fexpc = c.c;
out.statav5fexpr2 = gof.rsquare;
out.statav5fexpadjr2 = gof.adjrsquare;
out.statav5fexprmse = gof.rmse;

out.statav5diff = abs((stats(end,10)-stats(1,10))/stats(end,10));
out.statav5hp = gethp(stats(:,10));

% 11) swss5_1
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[-stats(end,11) -0.1 stats(end,11)]);
f = fittype('a*exp(b*x)+c','options',s);
[c, gof] = fit(r,stats(:,11),f);
out.swss5_1fexpa = c.a;
out.swss5_1fexpb = c.b;
out.swss5_1fexpc = c.c;
out.swss5_1fexpr2 = gof.rsquare;
out.swss5_1fexpadjr2 = gof.adjrsquare;
out.swss5_1fexprmse = gof.rmse;

out.swss5_1diff = abs((stats(end,11)-stats(1,11))/stats(end,11));
out.swss5_1hp = gethp(stats(:,11));

    function out = doyourcalcthing(y,y_rand)
        % Cross Correlation to original signal
        xc = xcorr(y,y_rand,1,'coeff');
        xcn1 = xc(1);
        xc1 = xc(3);
        
        % Norm of differences between original and randomized signals
        d1 = norm(y-y_rand) / length(y);
        
        % Autocorrelation
        ac1 = CO_autocorr(y_rand,1);
        ac2 = CO_autocorr(y_rand,2);
        ac3 = CO_autocorr(y_rand,3);
        ac4 = CO_autocorr(y_rand,4);
        
        % Entropies
        shen = EN_entropies(y_rand,'shannon');
        sampen = EN_sampenc(y_rand,2,0.2);
        
        % Stationarity
        statav5 = SY_StatAv(y_rand,5,'seg');
        swss5_1 = SY_slidwin(y_rand,'std','std',5,1);
        
        out = [xcn1, xc1, d1, ac1, ac2, ac3, ac4, shen, sampen, statav5, swss5_1];
    end

	function thehp = gethp(v)
		if v(end) > v(1)
			thehp = find(v > 0.5*(v(end)+v(1)),1,'first');
		else
			thehp = find(v < 0.5*(v(end)+v(1)),1,'first'); % last?
		end
	end

end