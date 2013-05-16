function out = DN_kssimp(x,pinkbeanie,ange)

[f xi] = ksdensity(x); % the smoothed empirical distribution

switch pinkbeanie
    case 'numpeaks' % number of peaks
        df = diff(f);
        ddf = diff(df);
        sdsp = ddf(sgnchange(df));
        out = length(find(sdsp<-0.0002)); % 'large enough' maxima
    case 'max'
        out = max(f);
    case 'entropy'
        out = -sum(f(f>0).*log(f(f>0))*(xi(2)-xi(1)));
    case 'numcross' % number of times crosses a given constant value
        if nargin<3
            disp('error, specify a number')
            out = NaN; return
        end
        places = sgnchange(f-ange);
        out = length(places);
    case 'area'
        if nargin<3
            disp('error, specify a number'), out=NaN; return
        end
        out=sum(f(f<ange).*(xi(2)-xi(1))); % integral under this portion
	case 'arclength'  % path length
		if nargin<3
            disp('error, specify a number'), out=NaN; return
        end
        m = mean(x);
		fd = abs(diff(f(xi>m-ange & xi<m+ange))); % The integrand in the path length formula
		out = sum(fd.*(xi(2)-xi(1)));
    case 'asym'
        m = mean(x);
        out1 = sum(f(xi>m).*(xi(2)-xi(1)));
        out2 = sum(f(xi<m).*(xi(2)-xi(1)));
        out = out1/out2;
    case 'plsym'
        m = mean(x);
        out1 = sum(abs(diff(f(xi<m))).*(xi(2)-xi(1)));
        out2 = sum(abs(diff(f(xi>m))).*(xi(2)-xi(1)));
        out = out1/out2;
end


end