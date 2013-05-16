function out = DN_kscomp(x,pinkbeanie,ange)
% Performs certain simple statistics on the difference between the 
% empirical distribution from data x, the standard distribution specified
% by pinkbeanie, and the statistic specified by ange
% Ben Fulcher 2009

%% PREPROCESSING
% fit the standard distribution
xtry = linspace(-40,40,400);
switch pinkbeanie
    case 'norm'
%         mu=0;sig=1; % should be (0,1), normalized; NO; we use raw ts now
        [a b] = normfit(x); ftry=normpdf(xtry,a,b);
    case 'ev'
        a = evfit(x); ftry=evpdf(xtry,a(1),a(2));
    case 'uni'
        [a b] = unifit(x); ftry=unifpdf(xtry,a,b);
    case 'beta'
        % clumsily scale to the range (0,1)
        x = (x-min(x)+0.01*std(x))/(max(x)-min(x)+0.02*std(x));
        a = betafit(x); ftry=betapdf(xtry,a(1),a(2));
    case 'rayleigh'
        if any(x<0),out=NaN; return
        else % valid domain to fit a Rayleigh distribution
            a=raylfit(x); ftry=raylpdf(xtry,a);
        end
    case 'exp'
        if any(x<0),out=NaN; return
        else a=expfit(x); ftry=exppdf(xtry,a);
        end
    case 'gamma'
        if any(x<0),out=NaN; return
        else a=gamfit(x); ftry=gampdf(xtry,a(1),a(2));
        end
    case 'gp'
        disp('forget about gp fits. Too difficult.');
        out = NaN; return
%         if any(x<0),out=NaN; return
%         else a=gpfit(x);
%         end
    case 'logn'
        if any(x<0), out = NaN; return
        else a = lognfit(x); ftry = lognpdf(xtry,a(1),a(2));
        end
    case 'wbl'
        if any(x<0),out=NaN; return
        else a=wblfit(x); ftry=wblpdf(xtry,a(1),a(2));
        end
end
xtmafit=xtry(sgnchange(ftry-1E-5)); % find where it hits 1E-5
xtmafit=[floor(xtmafit(1)*10)/10 ceil(xtmafit(end)*10)/10];

% determine smoothed empirical distribution
[f xi]=ksdensity(x);
xi=xi(f>1E-6); % only keep values greater than 1E-6
xi=[floor(xi(1)*10)/10 ceil(xi(end)*10)/10];

% find appropriate range [x1 x2] that incorporates the full range of both
x1=min([xtmafit(1) xi(1)]);
x2=max([xtmafit(2) xi(end)]);

%% Inefficient, but easier to just rerun both over the same range
xi=linspace(x1,x2,500);
f=ksdensity(x,xi); % the smoothed empirical distribution
switch pinkbeanie
    case 'norm'
        ffit=normpdf(xi,a,b);
    case 'ev'
        ffit=evpdf(xi,a(1),a(2));
    case 'uni'
        ffit=unifpdf(xi,a,b);
    case 'beta'
%         if x1<0,x1=1E-5;end
%         if x2>1,x2=1-1E-5;end
        ffit=betapdf(xi,a(1),a(2));
    case 'rayleigh'
%         if x1<0,x1=0;end
        ffit=raylpdf(xi,a);
    case 'exp'
%         if x1<0,x1=0;end
        ffit=exppdf(xi,a);
    case 'gamma'
%         if x1<0,x1=0;end
        ffit=gampdf(xi,a(1),a(2));
    case 'logn'
%         if x1<0,x1=0;end
        ffit=lognpdf(xi,a(1),a(2));
    case 'wbl'
%         if x1<0,x1=0;end
        ffit=wblpdf(xi,a(1),a(2));
end

% now the two cover the same range in x

%% CALCULATING
switch ange
    case 'adiff'
        % returns absolute area between the curves
        out=sum(abs(f-ffit)*(xi(2)-xi(1)));
    case 'peaksepy'
        % returns the seperation (in x) between the maxima of each distrn
        max1=max(f);
        max2=max(ffit);
        out=max2-max1;
    case 'peaksepx'
        % returns the seperation (in x) between the maxima of each distrn
        [max1 i1]=max(f);
        [max2 i2]=max(ffit);
        out=xi(i2)-xi(i1);
    case 'olapint'
        % returns the overlap integral between the two curves
        out=sum(f.*ffit*(xi(2)-xi(1)));
    case 'relent'
        % returns the relative entropy of the two distributions
        out=sum(f.*log(f./ffit)*(xi(2)-xi(1)));
end

end