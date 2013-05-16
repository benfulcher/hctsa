function out = OL_bentest(y,p,n)
% remove pth percent of highest and lowest outliers (i.e., 2*p percent
% removed altogether);
% returns the ratio of the mean (n==1) or std (n==2) before and after doing
% this
% Ben Fulcher

switch n
    case 1
        % look at mean (assume z-scored)
        out = mean(y(y>prctile(y,p) & y<prctile(y,100-p)));
    case 2
        % look at std
        out = std(y(y>prctile(y,p) & y<prctile(y,100-p)))/std(y); % [std(y) should be 1]
    otherwise
        disp('error in OL_bentest')
        return
end


end