function out = SD_makesurrogates(x,surrmethod,nsurrs,extrap)
% Ben Fulcher 27/1/2011
% Method described relatively clearly in Guarin Lopez et al. (arXiv, 2010)
% Used bits of aaft code that references (and presumably was obtained from)
% D. Kugiumtzis,
% ìSurrogate data test for nonlinearity including monotonic
% transformations,î Phys. Rev. E, vol. 62, no. 1, 2000.

% x should be input time series; surrogates should be made
% based on this.
% surrmethod will be used

% INPUTS
% number of surrogates to generate
if nargin < 3 || isempty(nsurrs)
    nsurrs = 1; % just create a single surrogate
end

% any extra parameters (some methods require)
if nargin < 4
    extrap = [];
end


N = length(x); % length of the time series
out = zeros(N,nsurrs); % each column is a new surrogate
tic

switch surrmethod
    case 'RP'
        % Random Phase Surrogates
        % Surrogates maintain linear correlations in the data, but any
        % nonlinear structure is destroyed by the phase randomization
        
        fprintf(1,'Constructing %u surrogates using the Random Phase Method\n',nsurrs)
        fprintf(1,'Linear correlations are maintained but nonlinear structure will be destroyed ' ...
                'by the phase randomization\n')
        
        % lost a datapoint if odd
        if rem(N,2) == 0
          n2 = N/2;
        else
          n2 = (N-1)/2;
        end
       
        for surri = 1:nsurrs
            % (*) Compute Fourier Transform of x => z
            z = fft(x,2*n2);
            
            % (*) Randomize Phases
            zMag = abs(z); % magnitude
            zPhase = angle(z); % phase
            
            randphase = 2*pi*rand(n2-1,1); % compute random phases
            
            % ensure phi(1)=0, and all others are in [0,2*pi]
            % (not quite sure what the zPhase(n2+1) is there for)...
            % negative phases to ensure complex conjugates -- IFT will be
            % real.
            newPhase = [0; randphase; zPhase(n2+1); -flipud(randphase)];
            
            
            % zNew is like z, but with randomized phases:
            zNew = [zMag(1:n2+1)', flipud(zMag(2:n2))']' .* exp(newPhase .* 1i);
            
            % Transform back into the time domain
            xNew = real(ifft(zNew,N));
            out(:,surri) = xNew;
        end
        
    case 'AAFT'
        fprintf(1,'Constructing %u surrogates using the Amplitude Adjusted Fourier '...
            'Transform (AAFT) Method\n',nsurrs)
        fprintf(1,'Linear correlations are maintained but nonlinear structure will be destroyed ' ...
                'by the phase randomization. Amplitude Distribution is approximately maintained\n')
        
        % Sort and rank order the data
        [xSorted, ix] = sort(x);
        [~, xRO] = sort(ix); % rank ordered permutation
        
        % lost a datapoint if odd
        if rem(N,2) == 0
          n2 = N/2;
        else
          n2 = (N-1)/2;
        end
        
        for surri = 1:nsurrs
            % Rand order white Gaussian-distributed noise
            nSort = sort(randn(N,1));
            y = nSort(xRO); % sorted Guassian white noise reordered as x

            % ------- Apply the RP method applied to y:
            % (*) Compute Fourier Transform of y => z
            z = fft(y,2*n2);
            
            % (*) Randomize Phases
            zMag = abs(z); % magnitude
            zPhase = angle(z); % phase
            
            randphase = 2*pi*rand(n2-1,1); % compute random phases
            
            % ensure phi(1)=0, and all others are in [0,2*pi]
            % (not quite sure what the zPhase(n2+1) is there for)...
            % negative phases to ensure complex conjugates -- IFT will be
            % real.
            newPhase = [0; randphase; zPhase(n2+1); -flipud(randphase)];
            
            
            % zNew is like z, but with randomized phases:
            zNew = [zMag(1:n2+1)', flipud(zMag(2:n2))']' .* exp(newPhase .* 1i);
            
            % Transform back into the time domain
            % phase-randomized version of random noise rank-ordered as x
            yRP = real(ifft(zNew,N));
            
            % --------- rank order x with respect to yRP
           [~, ixyRP] = sort(yRP);
           [~, yRO] = sort(ixyRP);
           out(:,surri) = xSorted(yRO);
        end        
        
    case 'TFT'
        fprintf(1,'Constructing %u surrogates using the Truncated Fourier '...
            'Transform (TFT) Method.\n',nsurrs)
        fprintf(1,'Low Frequency phases are preserved, and high frequency phases will be ' ...
                'randomized. A way of dealing with non-stationarity.\n')
        if isempty(extrap)
            fprintf(1,'You haven''t specified a cut-off frequency!! Setting N/8\n')
            fc = round(N/8);
        else
            fc = extrap; % extra input is the frequency cut-off
            if fc < 1
                fc = N*fc;
            end
        end
            
        % lost a datapoint if odd
        if rem(N,2) == 0
          n2 = N/2;
        else
          n2 = (N-1)/2;
        end
       
        for surri = 1:nsurrs
            % (*) Compute Fourier Transform of x => z
            z = fft(x,2*n2);
            
            % (*) Randomize Phases
            zMag = abs(z); % magnitude
            zPhase = angle(z); % phase
            
            randphase = pi*rand(n2-1,1); % compute random phases in (0,pi)
            randphase(1:fc) = zPhase(1:fc);
            
            % ensure phi(1)=0, and all others are in [0,2*pi]
            % (not quite sure what the zPhase(n2+1) is there for)...
            % negative phases to ensure complex conjugates -- IFT will be
            % real.
            newPhase = [0; randphase; zPhase(n2+1); -flipud(randphase)];
            
            
            % zNew is like z, but with randomized phases:
            zNew = [zMag(1:n2+1)' flipud(zMag(2:n2))']' .* exp(newPhase .* 1i);
            
            % Transform back into the time domain
            xNew = real(ifft(zNew,N));
            out(:,surri) = xNew;
        end
    otherwise
        error('Unknown surrogate generation method ''%s''',surrmethod)
end

% Cute farewell message
fprintf(1,'Generated %u %s surrogates in %s.\n',nsurrs,surrmethod,BF_thetime(toc,1))


end