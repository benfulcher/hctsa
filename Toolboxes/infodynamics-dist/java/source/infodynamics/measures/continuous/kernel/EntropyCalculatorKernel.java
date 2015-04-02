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

import infodynamics.measures.continuous.EntropyCalculator;

/**
 * <p>Computes the differential entropy of a given set of observations
 *  (implementing {@link EntropyCalculator}, using box-kernel estimation.
 *  For details on box-kernel estimation, see Kantz and Schreiber (below).</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link EntropyCalculator},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #EntropyCalculatorKernel()}.</li>
 *  <li>Further properties are available, see {@link #setProperty(String, String)};</li>
 *  <li>An additional {@link #initialise(double)} option;</li>
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
public class EntropyCalculatorKernel implements EntropyCalculator {

	protected KernelEstimatorUniVariate svke = null;
	/**
	 * Number of observations supplied
	 */
	protected int totalObservations = 0;
	/**
	 * Whether we're in debug mode
	 */
	protected boolean debug = false;
	/**
	 * The supplied observations
	 */
	protected double[] observations;
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
	public static final double DEFAULT_EPSILON = 0.25;
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
	public EntropyCalculatorKernel() {
		svke = new KernelEstimatorUniVariate();
		svke.setDebug(debug);
		svke.setNormalise(normalise);
	}

	public void initialise() {
		initialise(kernelWidth);
	}

	/**
	 * Initialise the calculator for (re-)use, with a specific kernel width,
	 * and existing (or default) values of other parameters.
	 * Clears an PDFs of previously supplied observations.
	 *
	 * @param kernelWidth if {@link #NORMALISE_PROP_NAME} property has
	 *  been set, then this kernel width corresponds to the number of
	 *  standard deviations from the mean (otherwise it is an absolute value)
	 */
	public void initialise(double kernelWidth) {
		this.kernelWidth = kernelWidth;
		svke.initialise(kernelWidth);
	}

	public void setObservations(double observations[]) {
		this.observations = observations;
		svke.setObservations(observations);
		totalObservations = observations.length;
	}
	
	public double computeAverageLocalOfObservations() {
		double entropy = 0.0;
		for (int t = 0; t < observations.length; t++) {
			double prob = svke.getProbability(observations[t]);
			double cont = Math.log(prob);
			entropy -= cont;
			if (debug) {
				System.out.println(t + ": p(" + observations[t] + ")= " +
						prob + " -> " + (cont/Math.log(2.0)) + " -> sum: " +
						(entropy/Math.log(2.0)));
			}
		}
		return entropy / (double) totalObservations / Math.log(2.0);
	}
	
	public void setDebug(boolean debug) {
		this.debug = debug;
		if (svke != null) {
			svke.setDebug(debug);
		}
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
			svke.setNormalise(normalise);
		} else {
			// No property was set
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

}
