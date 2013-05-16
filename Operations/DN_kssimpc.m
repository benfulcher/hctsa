function out = DN_kssimpc(x,pinkbeanie)
% REQUIRES RAW TIME SERIES AS INPUT: x

[f xi] = ksdensity(x); % the smoothed empirical distribution
[fz xiz] = ksdensity(zscore(x)); % smoothed zscored empirical distribution

switch pinkbeanie
    case 'numpeaks' % number of peaks
        df = diff(f); ddf=diff(df); % original
        sdsp = ddf(sgnchange(df));
        out1 = length(find(sdsp<-0.0002)); % 'large enough' maxima
        
        df = diff(fz); ddf = diff(df); % zscored
        sdsp = ddf(sgnchange(df));
        out2 = length(find(sdsp<-0.0002)); % 'large enough' maxima
        
        out = out2/out1; % shouldn't be meaningful
    case 'max'
        out1 = max(f);
        out2 = max(fz);
        out = out2/out1; % ratio of zscored to original maximum
    case 'entropy'
        out1 = -sum(f.*log(f)*(xi(2)-xi(1)));
        out2 = -sum(fz.*log(fz)*(xiz(2)-xiz(1)));
        out = out2/out1; % ratio of zscored to original entropy
        
%% DON'T USE THE REST; LEAVE THEM HERE FOR COMPLETENESS
%     case 'numcross' % number of times crosses a given constant value
%         if nargin<3
%             disp('error, specify a number')
%             out=NaN; return
%         end
%         places=sgnchange(f-ange);
%         out=length(places);
%     case 'area'
%         if nargin<3
%             disp('error, specify a number')
%             out=NaN; return
%         end
%         out=sum(f(f<ange).*(xi(2)-xi(1))); % integral under this portion
end


end