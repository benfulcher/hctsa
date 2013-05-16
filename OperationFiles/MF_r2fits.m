
function out=MF_r2fits(y,nbins,whoa,fox)

switch whoa
	case 'norm'% gaussian
		[a gof output]=fit(dnx',dny','gauss1');
	case 'exp'  % exponential fit
		[a gof output]=fit(dnx',dny','exp1');
	case 'pl' % power law
		[a gof output]=fit(dnx',dny','power1');
	%% SUPERCEDED BY MF_mtlbfit
% 	case 'sin2' % sin2
% 		try [a gof]=fit([1:length(y)]',y,'sin2');
% 	    catch bug
% 	        disp('error fitting sin2 model')
% 	        out=NaN; return
% 	    end
% 	case 'fourier1' % fourier1
% 		try [a gof]=fit([1:length(y)]',y,'fourier1');
% 	    catch bug
% 	        disp(['error fitting forier1 model (101) .. '])
% 	        out=NaN; return
% 	    end
end
switch fox
    case 'r2' % rsqured
        out=gof.rsquare;
    case 'adjr2' % degrees of freedom-adjusted rsqured
        out=gof.adjrsquare;
    case 'rmse' % root mean square error
        out=gof.rmse;
    case 'resAC1' % autocorrelation of residuals at lag 1
        out=CO_autocorr(output.residuals,1);
    case 'resAC2' % autocorrelation of residuals at lag 2
        out=CO_autocorr(output.residuals,2);
    case 'resruns' % runs test on residuals -- outputs p-value
        out=HT_hyptests(output.residuals,'runstest');
    case 'reslbq' % lbq test of randomness on residuals -- p-value
        out=HY_hyptests(output.residuals,'lbq');
end