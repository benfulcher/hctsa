figure
subplot(2,3,1);
echo on
clc

% This is an introduction to tstool
% tstool is a matlab toolbox for time series analysis
% 
% The basic object of interest in such a package is
% a "signal". A signal is a combination of raw data values
% and a descriptive part wich gives information about the
% raw data values, for example what physical unit they have
% and at which samplerate they were measured. Most information
% in the descriptive part is optional and not needed to do
% computations with that signal.

% When a signal is displayed at the matlab prompt, parts of
% the descriptive information is shown.

% Hit a key to continue 
pause; clc
% There are several different ways to create a signal object :

% The first methods is to use a matlab array as raw data for
% the new signal object. This gives you the possibility to use 
% results of other matlab scripts/functions/toolboxes in an easy
% way in the tstool package.

tmp = sin(0:0.1:79.9);
s = signal(tmp')

% Hit a key to continue
pause; clc
% Once you have a signal, you can use 
% command "view" to display the signal :

view(s,7);

% Hit a key to continue
pause; clc
% Okay, but this could be done in matlab without any additional packages. 

% Let's go to some more sophisticated stuff :
% Assume your time series was sampled with 8000 Hz and the data values
% are measured in Volts.

s = signal(tmp', 8000);
s = setyunit(s, unit('V'));
subplot(2,3,1);
view(s,7);

% Have a look at the plot now !
% Hit a key to continue
pause; clc
% Now its time to do some computations with signals.
% Let's beginn with the power spectrum of the signal.

subplot(2,3,2);
view(spec(s),7);

% What about the maximum of the recently calculated spectrum ?

max(spec(s))

% Hit a key to continue
pause; clc

% Let's view the auto correlation function (acf) of the clicks signal

s2 = acf(s);
subplot(2,3,3);
view(s2,7);


% Where is the first zero crossing ?

firstzero(s2)


% Hit a key to continue
pause; clc

% Another way to create a signal object is to load data from
% an ASCII, WAVE, AU or netCDF file

clicks = signal('clicks.wav', 'wav');

subplot(2,3,4);
view(clicks,7);

% Hit a key to continue
pause; clc
% For a spectrum in dB of the recently loaded time series use :

subplot(2,3,5)
view(db(spec(clicks)),7);

% By dragging a rectangular region in the plot (holding
% left mouse button down) you can zoom in.
% To zoom out, press the right mouse button several times.


% Hit a key to continue
pause; clc

% Spectrogramm

subplot(2,3,6);
view(swap(spec2(clicks)),7);

% Hit a key to continue
pause; clc




