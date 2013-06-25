function p = HT_disttests(x,ange,mydistn,nbins)
% Ben Fulcher

%% First fit the distribution
switch mydistn
    case 'norm'
        [a b] = normfit(x);
    case 'ev'
        a = evfit(x);
    case 'uni'
        [a b] = unifit(x);
    case 'beta'
        % clumsily scale to the range (0,1)
        x = (x-min(x)+0.01*std(x))/(max(x)-min(x)+0.02*std(x));
        a = betafit(x);        % then fit
    case 'rayleigh'
        if any(x<0),p=NaN; return
        else % valid domain to fit a Rayleigh distribution
            a=raylfit(x);
        end
    case 'exp'
        if any(x<0),p=NaN; return
        else a=expfit(x);
        end
    case 'gamma'
        if any(x<0),p=NaN; return
        else a=gamfit(x);
        end
    case 'gp'
        disp('forget about gp fits. Too difficult.');
%         if any(x<0),p=NaN; return
%         else a=gpfit(x);
%         end
    case 'logn'
        if any(x<=0),p=NaN; return
        else a=lognfit(x);
        end
    case 'wbl'
        if any(x<=0),p=NaN; return
        else a=wblfit(x);
        end
end

switch ange
    case 'chi2gof' % PERFORM A CHI2 GOODNESS OF FIT TEST
        switch mydistn
            case 'norm'
                mycdf = {@normcdf,a,b};
            case 'ev'
                mycdf = {@evcdf,a(1),a(2)};
            case 'uni'
                mycdf = {@unifcdf,a,b};
            case 'beta'
                mycdf = {@betapdf,a(1),a(2)};
            case 'rayleigh'
                mycdf = {@raylcdf,a};
            case 'exp'
                mycdf = {@expcdf,a};
            case 'gamma'
                mycdf = {@gamcdf,a(1),a(2)};
            case 'gp'
%                 mycdf = {@gpcdf,a(1),a(2),theta};
            case 'logn'
                mycdf = {@logncdf,a(1),a(2)};
            case 'wbl'
                mycdf = {@wblcdf,a(1),a(2)};
        end
        [h, p] = chi2gof(x,'cdf',mycdf,'nbins',nbins);
        
        
    case 'ks' % KOLMOGOROV-SMIRNOV TEST
        x1 = sort(unique(round(x*1E6)/1E6)); % unique values to 6 places
        if x1(1) > min(x), x1 = [min(x); x1]; end
        if x1(end) < max(x), x1 = [x1; max(x)]; end
        if size(x1,1) < size(x1,2);
            x1 = x1';
        end
        switch mydistn
            case 'norm'
                mycdf = [x1, normcdf(x1,a,b)];
            case 'ev'
                mycdf = [x1, evcdf(x1,a(1),a(2))];
            case 'uni'
                mycdf = [x1, unifcdf(x1,a,b)];
            case 'beta'
                mycdf = [x1, betacdf(x1,a(1),a(2))];
            case 'rayleigh'
                mycdf = [x1, raylcdf(x1,a)];
            case 'exp'
                mycdf = [x1, expcdf(x1,a)];
            case 'gamma'
                mycdf = [x1, gamcdf(x1,a(1),a(2))];
%             case 'gp'
%                 mycdf = [x1, gpcdf(x1,a(1),a(2))];
            case 'logn'
                mycdf = [x1, logncdf(x1,a(1),a(2))];
            case 'wbl'
                mycdf = [x1, wblcdf(x1,a(1),a(2))];
        end
%         rnondec = find(diff(mycdf(:,2))<=0);
%         if ~isempty(rnondec) % cdf is not strictly decreasing
%             mycdf(rnondec,:) = []; % remove these entries
%         end
        
        [h, p] = kstest(x,mycdf);
        
    case 'lillie' % LILLIEFORS TEST
        if any(ismember({'norm','ev'},mydistn))
            [h, p] = lillietest(x,0.05,mydistn);
        elseif strcmp('exp',mydistn)
            if any(x < 0)
                p = NaN; return
            else
                [h, p] = lillietest(x,0.05,mydistn);
            end
        else
           p = NaN;
           fprintf(1,'***RETURNED AN UNEXPECTED NAN FOR LILLIEFORS TEST\n')
        end
        
end

end