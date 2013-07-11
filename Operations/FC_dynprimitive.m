function out = FC_dynprimitive(y,fmeth)
% Looping over the length of data to use for prediction using
% FC_primitive, returns a bunch of statistics at each iteration
% Ben Fulcher, 2009

stats_st = zeros(10,7);

switch fmeth
    case 'mean'
        ngr = (1:10)';
    case 'median'
        ngr = (1:2:19)';
end

for i = 1:length(ngr)
    switch fmeth
        case 'mean'
            outtmp = FC_primitive(y,'mean',ngr(i));
        case 'median'
            outtmp = FC_primitive(y,'median',ngr(i));
            % median needs more tweaking
    end
%     stats_st(i,1)=outtmp.meanerr;
    stats_st(i,1) = outtmp.rmserr;
%     stats_st(i,2)=outtmp.meanabserr;
    stats_st(i,2) = outtmp.stderr;
    stats_st(i,3) = outtmp.sws;
    stats_st(i,4) = outtmp.swm;
%     stats_st(i,7)=outtmp.gofnadjr2;
    stats_st(i,5) = outtmp.ac1;
    stats_st(i,6) = outtmp.ac2;
%     stats_st(i,9)=outtmp.taures;
%     stats_st(i,10)=outtmp.tauresrat;
end

% plot(stats_st)
% out=stats_st;


% Get statistics from the shapes of the curves

% (1) root mean square error
% (i) (expect error to decrease with increasing forecast window?:)
out.rmserr_chn = mean(diff(stats_st(:,1)))/(range(stats_st(:,1)));
out.rmserr_meansgndiff = mean(sign(diff(stats_st(:,1))));

% (ii) is there a peak?
if out.rmserr_chn < 1; % on the whole decreasing, as expected
    wigv = max(stats_st(:,1));
    wig = find(stats_st(:,1) == wigv,1,'first');
    if wig~=1 && stats_st(wig-1,1)>wigv
        wig = NaN; % maximum is not a local maximum; previous value exceeds it
    elseif wig~=length(ngr) && stats_st(wig+1,1)>wigv
        wig = NaN; % maximum is not a local maximum; the next value exceeds it
    end
else
    wigv = min(stats_st(:,1));
    wig = find(stats_st(:,1) == wigv,1,'first');
    
    if wig~=1 && stats_st(wig-1,1)<wigv
        wig = NaN; % maximum is not a local maximum; previous value exceeds it
    elseif wig~=length(ngr) && stats_st(wig+1,1)<wigv
        wig = NaN; % maximum is not a local maximum; the next value exceeds it
    end
end
if ~isnan(wig)
    out.rmserr_peakpos = wig;
    out.rmserr_peaksize = wigv/mean(stats_st(:,1));
else % put NaNs in all the outputs
    out.rmserr_peakpos = NaN;
    out.rmserr_peaksize = NaN;
end

% (2) std of error
% this should correlate with (1), so let's just do a cross correlation
barrow = xcorr(stats_st(:,1),stats_st(:,2),1,'coeff');
out.xcorrstdrmserr = barrow(end); % xcorr at lag 1

% (3) Sliding Window Stationarity
out.sws_chn = mean(diff(stats_st(:,3)))/(range(stats_st(:,3)));
out.sws_meansgndiff = mean(sign(diff(stats_st(:,3))));
out.sws_stdn = std(stats_st(:,3))/range(stats_st(:,3));

% fit exponential decay
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[range(stats_st(:,3)) -0.5 min(stats_st(:,3))]);
f = fittype('a*exp(b*x)+c','options',s);
[c, gof] = fit(ngr,stats_st(:,3),f);
out.sws_fexp_a = c.a;
out.sws_fexp_b = c.b; % this is important
out.sws_fexp_c = c.c;
out.sws_fexp_r2 = gof.rsquare; % this is more important!
out.sws_fexp_adjr2 = gof.adjrsquare;
out.sws_fexp_rmse = gof.rmse;

% (4) sliding window mean
out.swm_chn = mean(diff(stats_st(:,4)))/(range(stats_st(:,4)));
out.swm_meansgndiff = mean(sign(diff(stats_st(:,4))));
out.swm_stdn = std(stats_st(:,4))/range(stats_st(:,4));

% (5) AC1
out.ac1_chn = mean(diff(stats_st(:,5)))/(range(stats_st(:,5)));
out.ac1_meansgndiff = mean(sign(diff(stats_st(:,5))));
out.ac1_stdn = std(stats_st(:,5))/range(stats_st(:,5));

% (6) AC2
out.ac2_chn=mean(diff(stats_st(:,6)))/(range(stats_st(:,6)));
out.ac2_meansgndiff=mean(sign(diff(stats_st(:,6))));
out.ac2_stdn=std(stats_st(:,6))/range(stats_st(:,6));

end