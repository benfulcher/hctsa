function Xs = DVV_surrogate(X, Ns)
% Generates surrogate data for real and complex signals
% 
% Surrogate data generation using iterated amplitude adjusted fourier
% transform (iAAFT) method for real-valued series, and complex iAAFT method
% for complex series.
%
%
% USAGE:    Xs = surrogate(X, Ns)
%
% INPUTS:
% X:        Input time series (can be real-valued or complex)
% Ns:       No. of surrogates to be generated
%
% OUTPUTS:
% Xs:   Generated surrogates
%
%
%   A Delay Vector Variance (DVV) toolbox for MATLAB
%   (c) Copyright Danilo P. Mandic 2008
%   http://www.commsp.ee.ic.ac.uk/~mandic/dvv.htm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default Parameters
if (nargin<1)
    error('Not enough Input arguments');
end
if (nargin<2)
    Ns = 25;
end

% Initial conditions and parameter initializations
iter = 0;
max_it = 1000;               % Maximum iterations
error_threshold = 1e-5;
MSE_start = 100;
MSE = 1000;

% Makes input vector X a column vector
if (size(X,2) > size(X,1))
    X = X';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements univariate iAAFT algorithm for real-valued signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isreal(X))

    Xs = zeros(length(X),Ns);

    % Stores the Amplitude vector of the original time series
    X_amp = abs(fft(X));

    % Stores the sorted version(ascending) of the original time series
    X_sorted = sort(X);

    % Iteration process for multiple surrogate time series generation
    for a = 1:Ns

        % Random permutation of the original time series
        temp = randperm(length(X));
        X_random = X(temp);

        % Initializations
        r_prev = X_random;
        MSE_prev = MSE_start;

        % Iterate untill convergence condition is met or max iterations are reached
        while (abs(MSE-MSE_prev) > error_threshold && iter < max_it)

            MSE_prev = MSE;

            % Amplitude spectrum matching
            ang_r_prev = angle(fft(r_prev));
            s = ifft(X_amp .* exp(ang_r_prev.*sqrt(-1)));

            % Rank ordering in order to scale to original signal distribution
            [s_sort, Ind] = sort(s);
            r(Ind,:) = X_sorted;

            % Convergence Metric calculation
            MSE = mean(abs(X_amp - abs(fft(r))));

            r_prev = r;
            iter = iter+1;
        end

        %         iter
        Xs(:,a) = r;
        iter = 0;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements Complex iAAFT Algorithm for complex signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else

    % Initial conditions
    Xs = zeros(length(X),Ns);
    r_real = zeros(length(X),1);
    r_imag = zeros(length(X),1);
    r = zeros(length(X),1);
    r_prev = zeros(length(X),1);

    % Stores the Amplitude vector of the original time series
    X_amp = abs(fft(X));

    % Stores the sorted version(modulus, ascending), of the original time series
    C_sorted = sort(abs(X));

    % Sorted versions(real and imaginary), of original time series
    C_real_sorted = sort(real(X));
    C_imag_sorted = sort(imag(X));

    % Iteration process for multiple surrogate time series generation
    for a = 1:Ns

        % Random permutation of the original time series
        r_prev = X(randperm(length(X)));
        MSE_prev = MSE_start;

         % Iterate untill convergence condition is met or max iterations are reached
        while( abs(MSE-MSE_prev) > error_threshold && iter < max_it)

            MSE_prev = MSE;

            % Amplitude spectrum matching
            ang_r_prev = angle(fft(r_prev));
            s = ifft(X_amp .* exp(ang_r_prev.*sqrt(-1)));

            % Rank ordering real and imaginary parts of s to match that of original signal
            [temp , Ind_real] = sort(real(s));
            [temp , Ind_imag] = sort(imag(s));
            r_real(Ind_real) = C_real_sorted;
            r_imag(Ind_imag) = C_imag_sorted;

            r = complex(r_real, r_imag);

            % Rank ordering in order to scale to original signal distribution
            [temp, Ind_abs] = sort(abs(r));
            AUX1 = abs(temp);
            r(Ind_abs) = r(Ind_abs).*(abs(C_sorted)./(AUX1+(AUX1==0)));

            % Convergence Metric
            MSE = mean(abs(X_amp - abs(fft(r))));
            r_prev = r;
            iter = iter + 1;
        end

        %         iter
        Xs(:,a) = r;
        iter = 0;
    end

end
