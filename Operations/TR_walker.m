function out = TR_walker(y,wstyle,wparam)
% This function uses the idea of simulating an artificial walker that 
% somehow moves in response to the time series.
% Outputs compare the motion of the walker to the actual observed time series
% wstyle determines how the walker's trajectory is defined in relation to the
%        time series
% wparam sets the parameters for the run
% Ben Fulcher August 2009

N = length(y); % the length of the input time series, y


%% (1) Walk
w = zeros(N,1); % the walker's trajectory
switch wstyle
    case 'prop'
        % walker starts at zero and narrows the gap between its position
        % and the time series value at that point by the proportion given
        % in wparam, to give the value at the subsequent time step
        p = wparam;
        
        w(1) = 0; % start at zero
        for i = 2:N
            w(i) = w(i-1)+p*(y(i-1)-w(i-1));
        end
    case 'biasprop'
        % walker is biased in one or the other direction (i.e., prefers to
        % go up, or down). Requires a vector of inputs: [p_up, p_down]
        pup = wparam(1); pdown = wparam(2);
        
        w(1) = 0;
        for i = 2:N
            if y(i) > y(i-1) % time series increases
                w(i) = w(i-1)+pup*(y(i-1)-w(i-1));
            else
                w(i) = w(i-1)+pdown*(y(i-1)-w(i-1));
            end
        end
    case 'momentum'
        % walker moves as if it had inertia from the previous time step,
        % i.e., it 'wants' to move the same amount; the time series acts as
        % a force changing its motion
        m = wparam(1); % 'inertial mass'
%         F=wparam(2); % weight of 'force' from time series
        
        w(1) = y(1);
        w(2) = y(2);
        for i = 3:N
            w_inert = w(i-1)+(w(i-1)-w(i-2));
%             w(i)=w_inert+(y(i-1)-w(i-1))/m; % dissipative term
            w(i) = w_inert+(y(i)-w_inert)/m; % dissipative term
            % equation of motion (s-s_0=ut+F/m*t^2)
            % where the 'force' F is the change in the original time series
            % at that point
            
        end
        
    case 'runningvar'
        % walker moves with momentum defined by amplitude of past values in
        % a given length window
        m = wparam(1); % 'inertial mass'
        wl = wparam(2); % window length
        
        w(1) = y(1);
        w(2) = y(2);
        for i = 3:N
            w_inert = w(i-1) + (w(i-1)-w(i-2));
            w_mom = w_inert + (y(i)-w_inert)/m; % dissipative term from time series
            if i > wl
                w(i) = w_mom*(std(y(i-wl:i))/std(w(i-wl:i))); % adjust by local standard deviation
            else
                w(i) = w_mom;
            end
        end
    otherwise
        error('Unknown method ''%s'' for simulating walker on the time series', wstyle)
end

%% % PLOT WALKER AND ORIGINAL TIME SERIES TOGETHER:
% lw = 1;
% figure('color','w'); hold on;
% c = bengetcmap('set1',3,1);
% plot(y,'.-k','LineWidth',lw); % original time series
% plot(w,'.-','color',c{1},'LineWidth',lw); % walker
% plot([1,length(w)],ones(2,1)*mean(w),'color',c{2},'LineWidth',2); % mean
% % running variance:
% stds = ones(N,2)*NaN;
% for i = wl+1:N
%     stds(i,1) = std(y(i-wl:i));
%     stds(i,2) = std(w(i-wl:i));
% end
% % plot(stds(:,1),':r'); % this is the time series
% plot(stds(:,1)./stds(:,2),'color',c{3},'LineWidth',lw); % this is the adjustment factor
% % means = zeros(N,1);
% % for i = 1:N
% %     means(i) = mean(w(1:i));
% % end
% % plot(means,'g')
% % plot(y-w,'m'); % residual
% legend('y','walker','mean: walker','localvariancefactor','accumulative walker mean')


%% (2) Statistics on walk
% (i) The walk itself
out.w_mean = mean(w);
out.w_median = median(w);
out.w_std = std(w);
out.w_ac1 = CO_autocorr(w,1);
out.w_ac2 = CO_autocorr(w,2);
out.w_tau = CO_fzcac(w);
out.w_min = min(w);
out.w_max = max(w);
out.w_propzcross = sum(w(1:end-1).*w(2:end) < 0)/(N-1);
% fraction of time series length that walker crosses time series

% (ii) Differences between the walk at signal
out.sw_meanabsdiff = mean(abs(y-w));
out.sw_taudiff = CO_fzcac(y) - CO_fzcac(w);
out.sw_stdrat = std(w)/std(y);
out.sw_ac1rat = out.w_ac1/CO_autocorr(y,1);
out.sw_minrat = min(w)/min(y);
out.sw_maxrat = max(w)/max(y);
out.sw_propcross = sum((w(1:end-1)-y(1:end-1)).*(w(2:end)-y(2:end)) < 0)/(N-1);
% fraction of time series length that walker crosses time series

% test from same distribution: Ansari-Bradley test
% (this may not be valid given the dependence of w on y, and the
% properties of the null hypothesis itself... But this is the name of the game!)
[h, pval, stats] = ansaribradley(w,y);
out.sw_ansarib_pval = pval; % p-value from the test
% out.sw_ansarib_W = stats.W; % W (test statistic)
% out.sw_ansarib_Wstar = stats.Wstar; % Approximate normal statistic
% test statistics are length dependent. Remove.

r = linspace(min(min(y),min(w)),max(max(y),max(w)),200); % make range of ksdensity uniform across all subsegments
dy = ksdensity(y,r); dw = ksdensity(w,r); % the kernel-smoothed distributions
out.sw_distdiff = sum(abs(dy-dw));

% (iii) Looking at residuals
res = w-y;
[h, pval] = runstest(res); % runs test
out.res_runstest = pval;
out.res_swss5_1 = SY_slidwin(res,'std','std',5,1); % sliding window stationarity
out.res_ac1 = CO_autocorr(res,1); % auto correlation at lag-1

end