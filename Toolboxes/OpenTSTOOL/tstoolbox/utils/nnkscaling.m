function se = nnkscaling(params)

% Mean square error function of scaling of neighbor distances
% ref. van de WATER and Schram : Generalized dimensions from near-neighbor information
% Phys Rev A Vol 37 Nr. 8 1988

global smallgamma k_values dl


try
	a = params(1);
	D_ = params(2);

	sc = a * (D_ * exp(gammaln(k_values+smallgamma/D_)-gammaln(k_values))).^(1/smallgamma);

	v = sc - dl;		% difference against measured values

	se = mean(v.*v);
catch
	se = 1e10 * mean(dl .* dl);
end

