function out = CO_acf_fit(y,whfn)
% Analyzes the decay properties of the autocorrlation function.
% The type of autocorrelation estimate is specified by whfn:
% whfn = 'acf' (autocorrelation function), 'glscf' (generalized linear self-correlation function)
% Keeps calculating until crosses zero, returns this lag (maximum = 400)
% Ben Fulcher

N = length(y);

switch whfn
	%% ACF
	case 'acf'
		% Calculate the function up to 12*tau (but not as much as N/2)
		tau = CO_fzcac(y);
		r = 0:1:min(12*tau,N/2);
		nr = length(r);
		acf = zeros(nr,1);
		for i = 1:nr;
			acf(i) = CO_autocorr(y,r(i));
		end

		% Statistics on Extrema		
		% Get the extrema
		rextrema = [0; find(diff(acf(1:end-1)).*diff(acf(2:end))<0)] + 1;
		extrema = acf(rextrema);
% 		hold off;plot(r,acf); hold on; plot(r(rextrema),extrema,'or');
		
		% Timmers calculated difference between absolute values of first two extrema of the acf:
		if length(extrema)>1;
			out.diff_ext12 = abs(extrema(2)-extrema(1));
			out.diff_ext12_tau = r(rextrema(2)) - r(rextrema(1));
		else out.diff_ext12 = NaN; out.diff_ext12_tau = NaN;
		end
		out.mean_ext_diff = mean(abs(diff(extrema)));
		out.nextrema = length(extrema);
	
		% Fit linear function to abs(extrema)
		p = polyfit(r(rextrema),extrema',1);
		out.linfit_ext_p1 = p(1);
		out.linfit_ext_p2 = p(2);
		extpred = polyval(p,r(rextrema));
		out.meanabserr = mean(abs(extpred-extrema'));
		
        %% Fits to the acf
		% Fit sinusoid across this range
		[cfun,gof] = fit(r,y,'sin1'); % fits form: a1*sin(b1*x+c1)
        keyboard
        out.sin1_gof_r2 = gof.rsquare;
        out.sin1_a1 = cfun.a1;
        out.sin1_b1 = cfun.b1;
		
        % Fit exponential decay to acf across this range
        s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 -1]);
        f = fittype('a*exp(b*x)','options',s);
        [c,gof] = fit(r(1:tau),y(1:tau),f);
        
        out.fitexpa = c.a;
        out.fitexpb = c.b;
        out.fitexpr2 = gof.rsquare;
        out.fitexpadjr2 = gof.adjrsquare;
        out.fitexprmse = gof.rmse;
        
		keyboard
	
	case 'glscf'
		glscfs = zeros(400,1);

		for i = 1:400
			tau = i;
			y1 = abs(y(1:end-tau));
			y2 = abs(y(1+tau:end));

			glscfs(i) = (mean((y1.^alpha).*(y2.^beta)) - mean(y1.^alpha)*mean(y2.^beta)) / ...
			 		( sqrt(mean(y1.^(2*alpha)) - mean(y1.^alpha)^2) * sqrt(mean(y2.^(2*beta)) - mean(y2.^beta)^2) );

			if i > 1 && glscfs(i)*glscfs(i-1)<0,
				% draw a straight line between these two and look at where hits zero
				out = i-1 + glscfs(i)/(glscfs(i)-glscfs(i-1));
				return;
			end
		end
		out = i;
end


end