figure
subplot(2,3,1);
echo on
clc

s = signal('colpitts.dat','ascii')

% Hit a key to continue
pause; clc
% colpitts is a signal from an electronical oscillator that shows
% nonlinear deterministic behaviour

view(s,7);

% Hit a key to continue
pause; clc

% Lets find a good choice for a delay-time
% by using the first minimum of the auto mutual information function
a = amutual(s,32);

subplot(2,3,2);
view(a,7);

% Hit a key to continue
pause; clc

% Now we need to know the minimal embedding dimension for the colpitts signal

c = cao(s,8,4,3,1000);
subplot(2,3,3);
view(c,7);

% Hit a key to continue
pause; clc

% Now do a time-delay reconstruction of the colpitts signal
e = embed(s, 3, 4);

subplot(2,3,4);
view(embed(cut(s,1,1000,1600),3,4), 7);

% Hit a key to continue
pause; clc
% What's the correlation dimension of the reconstructed data set ?

subplot(2,3,5)
view(corrsum2(e, [2000 100 5000], 0.1, 100, 32),7);

% Hit a key to continue
pause; clc

% And what about it's largest lyapunov exponent

subplot(2,3,6);
view(largelyap(e, 1000, 300, 40, 2),7);

% Hit a key to continue
pause; clc




