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

package infodynamics.measures.mixed.gaussian;

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.gaussian.EntropyCalculatorMultiVariateGaussian;
import infodynamics.measures.mixed.MutualInfoCalculatorMultiVariateWithDiscrete;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

/**
 * <p>Computes the differential mutual information between a given multivariate set of
 *  observations
 *  (<i>assuming that the probability distribution function for these observations is
 *  a multivariate Gaussian distribution</i>) 
 *  and a discrete variable.</p>
 *  
 * <p>This is done by examining the conditional probability distribution for
 *  the multivariate continuous variable C (given the discrete
 *  variable D) against the probability distribution for C:
 *  MI(C;D) := H(C) - H(C|D).</p>
 *  
 * <p><b>CAVEAT EMPTOR</b>: The real question this type of calculation asks
 * is to what extent does knowing the value of the discrete variable reduce
 * variance in the continuous variable(s). Indeed, it does not seem to behave
 * as we would normally expect a mutual information calculation: if we add more
 * continuous variables in, it increases in spite of redundancy between these
 * variables. TODO Further exploration should take place here ...</p>
 * 
 * <p>These calculators are <b>EXPERIMENTAL</b> -- not properly tested,
 * and not well documented. The intended calling pattern is similar to
 * {@link ConditionalMutualInfoCalculatorMultiVariate}
 * </p>
 * 
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct {@link #MutualInfoCalculatorMultiVariateWithDiscreteGaussian()}</li>
 *		<li>{@link #initialise(int, int)}</li>
 * 		<li>Set properties using {@link #setProperty(String, String)}</li>
 * 		<li>Provide the observations to the calculator using:
 * 			{@link #setObservations(double[][], int[])}, or
 * 			{@link #setCovariances(double[][], double[][])}.</li>
 * 		<li>Compute the required information-theoretic results, primarily:
 * 			{@link #computeAverageLocalOfObservations()} to return the average differential
 *          entropy based on either the set variance or the variance of
 *          the supplied observations; or other calls to compute
 *          local values or statistical significance.</li>
 * 	</ol>
 * </p>
 * 
 * @see Differential entropy for Gaussian random variables defined at 
 *      {@link http://mathworld.wolfram.com/DifferentialEntropy.html}
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class MutualInfoCalculatorMultiVariateWithDiscreteGaussian implements
		MutualInfoCalculatorMultiVariateWithDiscrete, Cloneable {

	/**
	 * Entropy calculator applied to the whole set of continuous data
	 */
	protected EntropyCalculatorMultiVariateGaussian entCalc;
	/**
	 * Entropy calculators applied to the set of continuous data
	 * associated with each discrete value
	 */
	protected EntropyCalculatorMultiVariateGaussian[] entCalcForEachDiscrete;
	
	/**
	 * Keep a copy of the discrete observations, to enable us to compute the
	 * statistical significance later.
	 */
	protected int[] discreteObservations;
	
	/**
	 * Keep a copy of the continuous observations, to enable us to compute the
	 * statistical significance later.
	 */
	protected double[][] continuousObservations;
	
	/**
	 * Number of supplied observations
	 */
	protected int totalObservations = 0;

	/**
	 * The number of possible discrete states
	 */
	protected int base = 0;

	/**
	 * Whether to print extra debug messages
	 */
	protected boolean debug = false;
	
	/**
	 * The last computed average MI value
	 */
	protected double lastAverage = 0;
	
	public MutualInfoCalculatorMultiVariateWithDiscreteGaussian() {
		entCalc = new EntropyCalculatorMultiVariateGaussian();
		entCalcForEachDiscrete = null;
	}
	
	public void initialise(int dimensions, int base) throws Exception {
		totalObservations = 0;
		lastAverage = 0;
		discreteObservations = null;
		this.base = base;
		entCalc.initialise(dimensions);
		entCalcForEachDiscrete = new EntropyCalculatorMultiVariateGaussian[base];
		for (int b = 0; b < base; b++) {
			entCalcForEachDiscrete[b] = new EntropyCalculatorMultiVariateGaussian();
			// If any properties relevant for these calculators were set in
			//  setProperty then we should set them here
			entCalcForEachDiscrete[b].initialise(dimensions);
		}
	}

	/**
	 * <p>Set the required property to the given value.</p>
	 * 
	 * <p>At this stage, there are no settable properties for this calculator.</p>
	 * 
	 * @param propertyName name of property
	 * @param propertyValue value of property
	 */
	public void setProperty(String propertyName, String propertyValue) {
		// No properties for this calculator
	}

	public void setObservations(double[][] continuousObservations,
			int[] discreteObservations) throws Exception {
		if (continuousObservations.length != discreteObservations.length) {
			throw new Exception("Observations are not of the same length");
		}
		// Set the complete set of observations:
		//  (this will pick up any errors in the dimensions of the continuous
		//   observations)
		entCalc.setObservations(continuousObservations);
		this.continuousObservations = continuousObservations;
		// Set the observations corresponding to each discrete value:
		setDiscreteData(continuousObservations, discreteObservations);
		totalObservations = continuousObservations.length;
	}

	protected void setDiscreteData(double continuousObservations[][],
			int discreteObservations[]) throws Exception {
		int totalNumberOfSuppliedObservations = 0;
		for (int b = 0; b < base; b++) {
			// Extract the observations for when this base value occurs:
			double[][] obsForThisDiscValue = MatrixUtils.extractSelectedPointsMatchingCondition(
					continuousObservations, discreteObservations, b);
			// Set the observations for each discrete value:
			entCalcForEachDiscrete[b].setObservations(obsForThisDiscValue);
			totalNumberOfSuppliedObservations += obsForThisDiscValue.length;
		}
		// Check that all of the supplied observations were extracted corresponding
		//  to one of the allowed discrete values
		if (totalNumberOfSuppliedObservations != discreteObservations.length) {
			throw new Exception("Some values in discreteObservations were not in the range 0..base-1");
		}
		this.discreteObservations = discreteObservations;
	}

	/**
	 * Compute the average mutual information from the previously supplied
	 *  observations via {@link #setObservations(double[][], int[])}
	 * 
	 * @return a scalar for the average mutual information in nats (NOT bits)
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations() throws Exception {
		// The average mutual information can be expressed
		//  as a difference between the entropy of the continuous observations
		//  and the conditional entropy of the continuous given the
		//  discrete observations:
		//		I(C;D) = H(C) - H(C|D)
		// Subtract that conditional entropy from the
		//  entropy of all observations:
		lastAverage = entCalc.computeAverageLocalOfObservations()
				- computeAverageLocalConditionalEntropyOfObservations();
		return lastAverage;
	}
	
	/**
	 * Compute the average conditional entropy of the continuous data given
	 * the discrete data (averaged over all discrete values)
	 * 
	 * @return average conditional entropy
	 * @throws Exception 
	 */
	protected double computeAverageLocalConditionalEntropyOfObservations() throws Exception {
		double meanConditionalEntropy = 0;
		for (int b = 0; b < base; b++) {
			double pOfB = (double) entCalcForEachDiscrete[b].getNumObservations() /
							(double) totalObservations;
			meanConditionalEntropy += pOfB *
					entCalcForEachDiscrete[b].computeAverageLocalOfObservations();
		}
		return meanConditionalEntropy;
	}

	/**
	 * Compute the local mutual information for each pair of continuous
	 *  and discrete values supplied here, using probability distribution
	 *  functions generated using the observations previously supplied
	 *  to the method {@link #setObservations(double[][], int[])}.
	 * 
	 * @param contStates observations of the joint continuous variables,
	 *  as specified in {@link #setObservations(double[][], int[])}
	 * @param discreteStates observations of the discrete variable,
	 *  as specified in {@link #setObservations(double[][], int[])}
	 * @return a time series of the local mutual information values
	 *  for each supplied observation, in nats (NOT bits)
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(
			double[][] contStates, int[] discreteStates) throws Exception {

		// The local mutual information can be expressed
		//  as a difference between the local 
		//  entropy of the continuous observations
		//  and the local conditional entropy of the continuous given the
		//  discrete observations:
		//		i(C;D) = h(C) - h(C|D)
		// First compute the local entropy for the continuous
		//  observations:
		double[] localValues = entCalc.computeLocalUsingPreviousObservations(contStates);
		// Next compute the local conditional entropes from each 
		//  conditional entropy calculator:
		double[][] localConditionalEntropies = new double[base][];
		for (int b = 0; b < base; b++) {
			// Extract the observations for when this base value occurs:
			double[][] obsForThisDiscValue = MatrixUtils.extractSelectedPointsMatchingCondition(
					contStates, discreteStates, b);
			//  and compute the local conditional entropies for these time points:
			localConditionalEntropies[b] =
					entCalcForEachDiscrete[b].
						computeLocalUsingPreviousObservations(obsForThisDiscValue);
		}
		// Now subtract the correct local conditional entropy from the local
		//  entropy of the continuous observations at the correct time points
		int[] nextTForDiscreteState = new int[base];
		for (int t = 0; t < contStates.length; t++) {
			// Find the local conditional entropy value at this time point:
			//  1. Grab the current discrete value
			int b = discreteStates[t];
			//  2. Pick out the next index in the local conditional entropies
			//    for this discrete value, nextTForDiscreteState[b], and pull
			//    out the local value at this index:
			double localCondEntropy =
					localConditionalEntropies[b][nextTForDiscreteState[b]];
			//  3. Update the next index for this discrete value:
			nextTForDiscreteState[b]++;
			// Now finalise the local MI at this point:
			localValues[t] -= localCondEntropy;
		}
		
		return localValues;
	}

	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		if (totalObservations == 0) {
			throw new Exception("Must have set observations before computing significance");
		}
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(totalObservations, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}

	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception {
		int numPermutationsToCheck = newOrderings.length;
		if (lastAverage == 0) {
			// This may execute even if the average was computed and found
			//  to be == 0; this won't hurt, just costs a tiny bit of execution
			//  time.
			computeAverageLocalOfObservations();
		}
		// Take a clone of the object to compute the MI of the surrogates:
		// (this is a shallow copy, it doesn't make new copies of all
		//  the arrays)
		MutualInfoCalculatorMultiVariateWithDiscreteGaussian miSurrogateCalculator =
				(MutualInfoCalculatorMultiVariateWithDiscreteGaussian) clone();
		
		double[] surrogateMeasurements = new double[numPermutationsToCheck];
		
		// Now compute the MI for each set of shuffled data:
		for (int i = 0; i < numPermutationsToCheck; i++) {
			// Generate a new re-ordered discrete data
			int[] shuffledDiscreteData = 
					MatrixUtils.extractSelectedTimePoints(
							discreteObservations, newOrderings[i]);
			// Re-initialise the data in the conditional entropy calculators:
			//  (in theory, we should call initialise and properly call
			//  setObservations(), but we know we can short circuit that from
			//  inside this calculator, and avoid recomputing covariances etc
			//  on the continuous data set)
			miSurrogateCalculator.setDiscreteData(
					continuousObservations, shuffledDiscreteData);
			// Compute the MI:
			surrogateMeasurements[i] = miSurrogateCalculator.entCalc.computeAverageLocalOfObservations() -
					miSurrogateCalculator.computeAverageLocalConditionalEntropyOfObservations();
			if (debug){
				System.out.println("New MI was " + surrogateMeasurements[i]);
			}
		}

		return new EmpiricalMeasurementDistribution(surrogateMeasurements, lastAverage);
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public double getLastAverage() {
		return lastAverage;
	}

	public int getNumObservations() {
		return totalObservations;
	}

	/**
	 * Clone the object - note: while it does create new cloned instances of
	 *  the {@link EntropyCalculatorMultiVariateGaussian} objects,
	 *  I think these only
	 *  have shallow copies to the data.
	 * This is enough though to maintain the structure across 
	 *  various {@link #computeSignificance(int)} calls.
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	protected Object clone() throws CloneNotSupportedException {
		MutualInfoCalculatorMultiVariateWithDiscreteGaussian theClone = 
				(MutualInfoCalculatorMultiVariateWithDiscreteGaussian) super.clone();
		// Now assign clones of the EntropyCalculatorMultiVariateGaussian objects:
		theClone.entCalc =
				(EntropyCalculatorMultiVariateGaussian) entCalc.clone();
		if (entCalcForEachDiscrete != null) {
			theClone.entCalcForEachDiscrete = new EntropyCalculatorMultiVariateGaussian[base];
			for (int b = 0; b < base; b++) {
				theClone.entCalcForEachDiscrete[b] =
						(EntropyCalculatorMultiVariateGaussian)
							entCalcForEachDiscrete[b].clone();
			}
		}
		return theClone;
	}
}
