%tstoolbox/mex/chaosys
%   chaosys gives the user the possibility to compute time series data for
%   a couple of dynamical systems, among which are Lorenz, Chua, Roessler
%   etc. This routine is not meant as a replacement for Matlab's suite of
%   functions for solving ODEs, but as a fast way to generate some data
%   sets to evaluate the processing capabilities of TSTOOL. The
%   integration is done by an ODE solver using an Adams Pece scheme with
%   local extrapolation . It is at least faster than However, it is not
%   possible to extend the set of systems without recompiling chaosys.
%
%   Syntax:
%
%     * x = chaosys(length, stepwidth, initial_conditions, mode,
%       parameters)
%
%   Input arguments:
%
%     * length - number of samples to generate
%     * stepwidth - integration step size
%     * initial_conditions - vector of initial conditions
%     * mode:
%          + 0: Lorenz
%          + 1: Generalized Chua : Double Scroll
%          + 2: Generalized Chua : Five Scroll
%          + 3: Duffing
%          + 4: Roessler
%          + 5: Toda Oscillator
%          + 6: Van der Pol Oscillator
%          + 7: Pendulum
%       For an exact definition of the ODE systems, please refer to this
%       header file.
%     * parameters - vector of systems parameters. The order of the
%       parameters is exactly the same as in the constructors of the DGL
%       subclasses in the above file.
%
%   Output arguments:
%
%     * x contains the output of the integration, organized as matrix of
%       size samples by dim, where dim is the number of ODEs that define
%       the system
%
%   Example:
%
%x = chaosys(20000, 0.025, [0.1 -0.1 0.02], 0);
%plot(x(:,1));
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

