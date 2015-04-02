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

package infodynamics.measures.mixed;

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.measures.mixed.ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSource;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

/**
 * <p>This is an abstract calculator to
 *  compute the Conditional Mutual Information I(X;D|Z) 
 *  between a discrete variable D and a 
 *  vector of continuous variables X,
 *  conditioned on another vector of continuous variables Z.
 *  It implements common methods for child classes to build on</p>
 *  
 * <p>These calculators are <b>EXPERIMENTAL</b> -- not properly tested,
 * and not well documented. The intended calling pattern is similar to
 * {@link ConditionalMutualInfoCalculatorMultiVariate}
 * </p>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public abstract class ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSourceCommon
	implements ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSource {

	/**
	 * we compute distances to the kth neighbour
	 */
	protected double[][] continuousDataX;
	protected double[][] conditionedDataZ;
	protected int[] discreteData;
	protected int[] counts;
	protected int base;
	protected int dimensionsContinuous;
	protected int dimensionsConditional;
	protected double[] meansX;
	protected double[] stdsX;
	protected double[] meansZ;
	protected double[] stdsZ;
	protected boolean debug;
	protected double condMi;
	protected boolean miComputed;
	protected int totalObservations;
	
	public static final String PROP_NORMALISE = "NORMALISE";
	protected boolean normalise = true;

	public ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSourceCommon() {
	}

	/**
	 * Initialise the calculator.
	 * 
	 * @param dimensions number of joint continuous variables
	 * @param base number of discrete states
	 * @param dimensionsCond the number of joint continuous variables
	 *     to condition on
	 */
	public void initialise(int dimensions, int base, int dimensionsCond) {
		condMi = 0.0;
		miComputed = false;
		continuousDataX = null;
		discreteData = null;
		this.base = base;
		dimensionsContinuous = dimensions;
		dimensionsConditional = dimensionsCond;
	}

	/**
	 * Sets properties for the calculator.
	 * Valid properties include:
	 * <ul>
	 *  <li>{@link #PROP_NORMALISE} - whether to normalise the individual
	 *      variables (true by default)</li>
	 * </ul>
	 * 
	 * @param propertyName
	 * @param propertyValue
	 */
	public void setProperty(String propertyName, String propertyValue) {
		if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
			normalise = Boolean.parseBoolean(propertyValue);
		}
	}

	public void startAddObservations() {
		throw new RuntimeException("Not implemented yet");
	}

	public void finaliseAddObservations() throws Exception {
		meansX = MatrixUtils.means(continuousDataX);
		stdsX = MatrixUtils.stdDevs(continuousDataX, meansX);
		meansZ = MatrixUtils.means(conditionedDataZ);
		stdsZ = MatrixUtils.stdDevs(conditionedDataZ, meansZ);
		if (normalise) {
			// Take a copy since we're going to normalise it
			continuousDataX = MatrixUtils.normaliseIntoNewArray(continuousDataX);
			conditionedDataZ = MatrixUtils.normaliseIntoNewArray(conditionedDataZ);
		}
		// count the discrete states:
		counts = new int[base];
		for (int t = 0; t < discreteData.length; t++) {
			counts[discreteData[t]]++;
		}
		totalObservations = discreteData.length;
	}

	public void addObservations(double[][] continuousObservations,
			int[] discreteObservations, double[][] conditionedObservations)
			throws Exception {
		throw new RuntimeException("Not implemented yet");
	}
	
	public void setObservations(double[][] continuousObservations,
			int[] discreteObservations, double[][] conditionedObservations)
			throws Exception {
		
		if ((continuousObservations.length != discreteObservations.length) ||
			(continuousObservations.length != conditionedObservations.length)) {
			throw new Exception("Time steps for observations2 " +
					discreteObservations.length + " does not match the length " +
					"of observations1 " + continuousObservations.length +
					" and of conditionedObservations " + conditionedObservations.length);
		}
		if (continuousObservations[0].length == 0) {
			throw new Exception("Computing MI with a null set of data");
		}
		if (conditionedObservations[0].length == 0) {
			throw new Exception("Computing MI with a null set of conditioned data");
		}
		continuousDataX = continuousObservations;
		discreteData = discreteObservations;
		conditionedDataZ = conditionedObservations;
		
		finaliseAddObservations();
		
	}

	/**
	 * Compute the significance of the mutual information of the previously supplied observations.
	 * We destroy the p(x,d|z) correlations, while retaining the p(x|z) marginals, to check how
	 *  significant this conditional mutual information actually was.
	 *  
	 * This is in the spirit of Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128
	 *  which was performed for Transfer entropy.
	 * 
	 * @param reorderDiscreteVariable boolean for whether to reorder the discrete variable (true)
	 *    or the continuous variable (false)
	 * @param numPermutationsToCheck
	 * @return the proportion of conditional MI scores from the distribution which have higher or equal MIs to ours.
	 */
	public synchronized EmpiricalMeasurementDistribution computeSignificance(
			boolean reorderDiscreteVariable, int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		// Use continuousDataX length (all variables have same length) even though
		//  we may be randomising the other variable:
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(
				continuousDataX.length, numPermutationsToCheck);
		return computeSignificance(reorderDiscreteVariable, newOrderings);
	}
	
	/**
	 * <p>As per {@link #computeSignificance(boolean, int)} but supplies
	 *  the re-orderings of the observations of the named variable.</p>
	 * 
	 * <p>We provide a simple implementation which would be suitable for
	 *  any child class, though the child class may prefer to make its
	 *  own implementation to make class-specific optimisations.
	 *  Child classes must implement {@link java.lang.Cloneable}
	 *  for this method to be callable for them, and indeed implement
	 *  the clone() method in a way that protects their structure
	 *  from alteration by surrogate data being supplied to it.</p>
	 * 
	 * @param reorderDiscreteVariable boolean for whether to reorder the discrete variable (true)
	 *    or the continuous variable (false)
	 * @param newOrderings the specific new orderings to use
	 * @return the proportion of conditional MI scores from the distribution which have higher or equal MIs to ours.
	 */
	public EmpiricalMeasurementDistribution computeSignificance(
			boolean reorderDiscreteVariable, int[][] newOrderings) throws Exception {
		
		int numPermutationsToCheck = newOrderings.length;
		if (!miComputed) {
			computeAverageLocalOfObservations();
		}
		
		// Take a clone of the object to compute the MI of the surrogates:
		// (this is a shallow copy, it doesn't make new copies of all
		//  the arrays - child classes should override this)
		ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSource miSurrogateCalculator =
				(ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSource) this.clone();
		
		double[] surrogateMeasurements = new double[numPermutationsToCheck];
		
		// Now compute the MI for each set of shuffled data:
		for (int i = 0; i < numPermutationsToCheck; i++) {
			if (reorderDiscreteVariable) {
				// Reorder the discrete variable
				int[] shuffledData = 
					MatrixUtils.extractSelectedTimePoints(discreteData, newOrderings[i]);
				// Perform new initialisations
				miSurrogateCalculator.initialise(
						dimensionsContinuous, base, dimensionsConditional);
				// Set new observations
				miSurrogateCalculator.setObservations(continuousDataX,
						shuffledData, conditionedDataZ);
			} else {
				// Reorder the continuous variable
				double[][] shuffledData = 
						MatrixUtils.extractSelectedTimePointsReusingArrays(
								continuousDataX, newOrderings[i]);
				// Perform new initialisations
				miSurrogateCalculator.initialise(
						dimensionsContinuous, base, dimensionsConditional);
				// Set new observations
				miSurrogateCalculator.setObservations(shuffledData,
						discreteData, conditionedDataZ);
			}
			// Compute the MI
			surrogateMeasurements[i] = miSurrogateCalculator.computeAverageLocalOfObservations();
			if (debug){
				System.out.println("New MI was " + surrogateMeasurements[i]);
			}
		}
		
		return new EmpiricalMeasurementDistribution(surrogateMeasurements, condMi);
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public double getLastAverage() {
		return condMi;
	}
	
	public int getNumObservations() {
		return totalObservations;
	}
}
