/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2012, Joseph T. Lizier
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package infodynamics.measures.continuous.kernel;

import infodynamics.measures.continuous.EntropyCalculatorMultiVariate;

/**
 * <p>Computes the differential entropy of a given set of observations
 *  (implementing {@link EntropyCalculatorMultiVariate}, using box-kernel estimation.
 *  For details on box-kernel estimation, see Kantz and Schreiber (below).</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link EntropyCalculatorMultiVariate},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #EntropyCalculatorMultiVariateKernel()}.</li>
 *  <li>Further properties are available, see {@link #setProperty(String, String)};</li>
 *  <li>An additional {@link #initialise(int, double)} option;</li>
 * </ul>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>H. Kantz and T. Schreiber, "Nonlinear Time Series Analysis".
 *   Cambridge, MA: Cambridge University Press, 1997.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class EntropyCalculatorMultiVariateKernel implements EntropyCalculatorMultiVariate {

	private KernelEstimatorMultiVariate mvke = null;
	/**
	 * Number of observations supplied
	 */
	private int totalObservations = 0;
	// private int dimensions = 0;
	/**
	 * Whether we're in debug mode
	 */
	private boolean debug = false;
	/**
	 * The supplied observations
	 */
	private double[][] observations = null;
	/**
	 * Last computed average entropy
	 */
	private double lastEntropy;
	/**
	 * Whether we normalise the incoming observations to mean 0,
	 * standard deviation 1.
	 */
	private boolean normalise = true;
	/**
	 * Property for whether we normalise the incoming observations to mean 0,
	 * standard deviation 1.
	 */
	public static final String NORMALISE_PROP_NAME = "NORMALISE";
	
	/**
	 * Default value for kernel width
	 */
	private static final double DEFAULT_EPSILON = 0.25;
	/**
	 * Kernel width
	 */
	private double kernelWidth = DEFAULT_EPSILON;
	/**
	 * Property name for the kernel width
	 */
	public static final String KERNEL_WIDTH_PROP_NAME = "KERNEL_WIDTH";
	/**
	 * Legacy property name for the kernel width
	 */
	public static final String EPSILON_PROP_NAME = "EPSILON";

	/**
	 * Construct an instance
	 */
	public EntropyCalculatorMultiVariateKernel() {
		mvke = new KernelEstimatorMultiVariate();
		mvke.setDebug(debug);
		mvke.setNormalise(normalise);
		lastEntropy = 0.0;
	}

	public void initialise(int dimensions) {
		initialise(dimensions, kernelWidth);
	}

	/**
	 * Initialise the calculator for (re-)use, with a specific 
	 * number of joint variables, specific kernel width,
	 * and existing (or default) values of other parameters.
	 * Clears an PDFs of previously supplied observations.
	 *
	 * @param dimensions number of joint variables
	 * @param kernelWidth if {@link #NORMALISE_PROP_NAME} property has
	 *  been set, then this kernel width corresponds to the number of
	 *  standard deviations from the mean (otherwise it is an absolute value)
	 */
	public void initialise(int dimensions, double kernelWidth) {
		this.kernelWidth = kernelWidth;
		mvke.initialise(dimensions, kernelWidth);
		// this.dimensions = dimensions;
		lastEntropy = 0.0;
	}

	public void setObservations(double observations[][]) {
		mvke.setObservations(observations);
		totalObservations = observations.length;
		this.observations = observations;
	}

	public double computeAverageLocalOfObservations() {
		double entropy = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob = mvke.getProbability(observations[b]);
			double cont = Math.log(prob);
			entropy -= cont;
			if (debug) {
				System.out.println(b + ": " + prob + " -> " + (-cont/Math.log(2.0)) + " -> sum: " + (entropy/Math.log(2.0)));
			}
		}
		lastEntropy = entropy / (double) totalObservations / Math.log(2.0);
		return lastEntropy;
	}

	
	public double[] computeLocalOfPreviousObservations() {
		return computeLocalUsingPreviousObservations(observations);
	}

	public double[] computeLocalUsingPreviousObservations(double states[][]) {
		double entropy = 0.0;
		double[] localEntropy = new double[states.length];
		for (int b = 0; b < states.length; b++) {
			double prob = mvke.getProbability(states[b]);
			double cont = -Math.log(prob);
			localEntropy[b] = cont;
			entropy += cont;
			if (debug) {
				System.out.println(b + ": " + prob + " -> " + (cont/Math.log(2.0)) + " -> sum: " + (entropy/Math.log(2.0)));
			}
		}
		entropy /= (double) totalObservations / Math.log(2.0);
		return localEntropy;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
		mvke.setDebug(debug);
	}

	public double getLastAverage() {
		return lastEntropy;
	}
	
	/**
	 * <p>Set properties for the kernel entropy calculator.
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 * 		<li>{@link #KERNEL_WIDTH_PROP_NAME} (legacy value is {@link #EPSILON_PROP_NAME}) --
	 * 			kernel width to be used in the calculation. If {@link #normalise} is set,
	 * 		    then this is a number of standard deviations; otherwise it
	 * 			is an absolute value. Default is {@link #DEFAULT_KERNEL_WIDTH}.</li>
	 * 		<li>{@link #NORMALISE_PROP_NAME} -- whether to normalise the incoming variable values
	 * 			to mean 0, standard deviation 1, or not (default false). Sets {@link #normalise}.</li>
	 * </ul> 
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		boolean propertySet = true;

		// TODO If we implement a dynamic correlation exclusion property,
		//  then we will need to call getProbability(double, int) instead of
		//  just getProbability(double) above.
		
		if (propertyName.equalsIgnoreCase(KERNEL_WIDTH_PROP_NAME) ||
				propertyName.equalsIgnoreCase(EPSILON_PROP_NAME)) {
			kernelWidth = Double.parseDouble(propertyValue);
		} else if (propertyName.equalsIgnoreCase(NORMALISE_PROP_NAME)) {
			normalise = Boolean.parseBoolean(propertyValue);
			mvke.setNormalise(normalise);
		} else {
			// No property was set
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	public int getNumObservations() throws Exception {
		return totalObservations;
	}

}
