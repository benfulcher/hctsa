% Cross Recurrence Plot Toolbox.
% Version 5.17  (R28.15) 14-Nov-2012
%
%   ace            - Finds optimal transformation and maximal correlation.
%   adjust         - Adjusts two two-column vectors.
%   arfit          - AR parameter estimation via Yule-Walker method.
%   choosecolormap - GUI for choosing a colormap.
%   corrgram       - Calculate windowed cross correlation between two signals.
%   crp            - Cross recurrence plot/ recurrence plot of given data series.
%   crp_big        - Cross recurrence plot/ recurrence plot of given data series for long data series.
%   crp2           - Cross recurrence plot/ recurrence plot of given vectors and the line of synchronization.
%   crpclean       - Removes the CRP toolbox.
%   crqa           - Recurrence quantification analysis.
%   crqad          - Recurrence quantification analysis diagonalwise.
%   crqad          - Recurrence quantification analysis diagonalwise for long data series.
%   dl             - Distribution and measures of diagonal RP structures.
%   entropy        - Entropy of a distribution.
%   fnn            - Dimension estimation by false nearest neighbours.
%   french         - French's flag color map.
%   hist2          - Two-dimensional histogram.
%   histn          - Multi-dimensional histogram.
%   jrp            - Joint recurrence plot/ recurrence plot of given data series.
%   jrqa           - Joint recurrence quantification analysis.
%   mcf            - Plots maximal correlation function.
%   mgui           - Starts a GUI for data analysis programmes.
%   mi             - Multi-dimensional mutual information.
%   migram         - Calculate windowed mutual information between two signals.
%   normalize      - Normalizes data series.
%   phasespace     - Embedding of time series in a phase space.
%   phasesynchro   - Indicator of phase synchronisation by means of recurrences.
%   pss            - Computes phase space size.
%   recons         - Reconstruct a time series from a recurrence plot.
%   rpde           - Computes the recurrence time entropy.
%   rrspec         - Tau-recurrence rate spectrum.
%   rtspec         - Recurrence time spectrum.
%   taucrp         - Creates a close returns plot.
%   trackplot      - Estimates the line of synchronization in a CRP. 
%   trafo          - Transforms data to a desired distribution.
%   tt             - Distribution and measures of vertical RP structures.
%   twinsurr       - Creates twin surrogates for statistical tests.
%   waitbar        - Display wait bar (the improved Mathworks waitbar).
%   winplot        - Windowed plot.
%   xcf            - Computes and plots crosscorrelation.
%
%
% Copyright (c) 2008-2010
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 1998-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% This toolbox is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License 
%        gpl - click here to view the GNU Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Modified at 14-Nov-2012 14:07:46 by MAKEINSTALL
