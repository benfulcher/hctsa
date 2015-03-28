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

import infodynamics.measures.continuous.TransferEntropyCalculator;
import infodynamics.measures.continuous.TransferEntropyCommon;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;

import java.util.Iterator;

/**
 * <p>Computes the differential transfer entropy (TE) between two uniivariate
 *  <code>double[]</code> time-series of observations
 *  using box-kernel estimation.
 *  For details on box-kernel estimation, see Kantz and Schreiber (below).</p>
 *  
 * <p>TE was defined by Schreiber (below).
 *  This estimator is realised here by extending
 *  {@link TransferEntropyCommon}.</p>
 *  
 * <p>***** This implementation is <i>simple</i>, involving O(n^2) comparisons.
 * The main TE box-kernel class, {@link TransferEntropyCalculatorKernel} is significantly
 * more efficient.
 * This class incurs a larger memory cost than {@link TransferEntropyCalculatorKernelPlainIterators}
 * but has a significantly lower cost in time than it (validate this, it has been
 * a long time since I wrote that!).
 * </p>
 * 
 * <p><b>Since this is not our main TE box-kernel class, the javadocs are
 * not necessarily well-maintained either ...</b> TODO Update javadocs here</p>
 * 
 * <p>Usage is as per the paradigm outlined for {@link TransferEntropyCalculator},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #TransferEntropyCalculatorKernel()}.</li>
 *  <li>Further properties are available, see {@link #setProperty(String, String)};</li>
 *  <li>An additional {@link #initialise(int, double)} option;</li>
 *  <li>Additional utility methods for computing other information-theoretic values
 *  are available here (e.g. {@link #computeAverageLocalOfObservationsWithCorrection()}) which
 *  can be called after all observations are supplied.</li>
 *  </ul>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. Schreiber, <a href="http://dx.doi.org/10.1103/PhysRevLett.85.461">
 * "Measuring information transfer"</a>,
 *  Physical Review Letters 85 (2) pp.461-464, 2000.</li>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href="http://dx.doi.org/10.1103/PhysRevE.77.026110">
 *  "Local information transfer as a spatiotemporal filter for complex systems"</a>
 *  Physical Review E 77, 026110, 2008.</li>
 * </ul>
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class TransferEntropyCalculatorKernelPlain
	extends TransferEntropyCommon implements TransferEntropyCalculator {

	// Keep joint vectors so we don't need to regenerate them
	protected double[][] destNextPastSourceVectors;
	
	// Indices into array of count results
	private static final int NEXT_PAST_SOURCE = 0;
	private static final int PAST_SOURCE = 1;
	private static final int NEXT_PAST = 2;
	private static final int PAST = 3;	
	
	private boolean normalise = true;
	public static final String NORMALISE_PROP_NAME = "NORMALISE";
	
	private boolean dynCorrExcl = false;
	private int dynCorrExclTime = 100;
	public static final String DYN_CORR_EXCL_TIME_NAME = "DYN_CORR_EXCL";
	
	/**
	 * Default value for epsilon
	 */
	public static final double DEFAULT_EPSILON = 0.25;
	/**
	 * Kernel width
	 */
	private double epsilon = DEFAULT_EPSILON;
	private double[] epsilons;
	public static final String EPSILON_PROP_NAME = "EPSILON";
	
	/**
	 * Creates a new instance of the kernel-estimate style transfer entropy calculator
	 *
	 */
	public TransferEntropyCalculatorKernelPlain() {
		super();
	}

	/**
	 * Initialises the calculator with the existing value for epsilon
	 * 
	 * @param k history length
	 */
	public void initialise(int k) throws Exception {
		initialise(k, epsilon);
	}
	
	/**
	 * Initialises the calculator
	 * 
	 * @param k history length
	 * @param epsilon kernel width
	 */
	public void initialise(int k, double epsilon) throws Exception {
		this.epsilon = epsilon;
		super.initialise(k); // calls initialise()
	}
	
	/**
	 * Initialise using default or existing values for k and epsilon
	 */
	public void initialise() {
		destNextPastSourceVectors = null;
		epsilons = new double[k + 2];		
	}

	/**
	 * Set properties for the transfer entropy calculator.
	 * These can include:
	 * <ul>
	 * 		<li>EPSILON_PROP_NAME</li>
	 * 		<li>NORMALISE_PROP_NAME</li>
	 * 		<li>DYN_CORR_EXCL_TIME_NAME</li>
	 * </ul> 
	 * 
	 * @param propertyName
	 * @param propertyValue
	 * @throws Exception 
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		if (propertyName.equalsIgnoreCase(EPSILON_PROP_NAME)) {
			epsilon = Double.parseDouble(propertyValue);
		} else if (propertyName.equalsIgnoreCase(NORMALISE_PROP_NAME)) {
			normalise = Boolean.parseBoolean(propertyValue);
		} else if (propertyName.equalsIgnoreCase(DYN_CORR_EXCL_TIME_NAME)) {
			dynCorrExclTime = Integer.parseInt(propertyValue);
			dynCorrExcl = (dynCorrExclTime > 0);
		} else {
			// try the superclass:
			super.setProperty(propertyName, propertyValue);
		}

		if (debug) {
			System.out.println("Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	/**
	 * Flag that the observations are complete, probability distribution functions can now be built.
	 *
	 */
	public void finaliseAddObservations() {
		// First work out the size to allocate the joint vectors, and do the allocation:
		totalObservations = 0;
		for (double[] destination : vectorOfDestinationObservations) {
			totalObservations += destination.length - k;
		}
		destNextPastSourceVectors = new double[totalObservations][k + 2];
		
		// Construct the joint vectors from the given observations
		int startObservation = 0;
		Iterator<double[]> iterator = vectorOfDestinationObservations.iterator();
		for (double[] source : vectorOfSourceObservations) {
			double[] destination = iterator.next();
			double[][] currentDestNextPastSourceVectors = makeJointVectorForNextPastSource(destination, source);
			MatrixUtils.arrayCopy(currentDestNextPastSourceVectors, 0, 0,
					destNextPastSourceVectors, startObservation, 0, currentDestNextPastSourceVectors.length, k + 2);
			startObservation += destination.length - k;
		}
		
		if (normalise) {
			// Adjust epsilon for each dimension
			for (int c = 0; c < k + 2; c++) {
				double std = MatrixUtils.stdDev(destNextPastSourceVectors, c);
				epsilons[c] = epsilon * std;
			}
		} else {
			for (int c = 0; c < k + 2; c++) {
				epsilons[c] = epsilon;
			}
		}
		
		// Store whether there was more than one observation set:
		addedMoreThanOneObservationSet = vectorOfDestinationObservations.size() > 1;
		if (addedMoreThanOneObservationSet && dynCorrExcl) {
			// We have not properly implemented dynamic correlation exclusion for
			//  multiple observation sets, so throw an error
			throw new RuntimeException("Addition of multiple observation sets is not currently " +
					"supported with property DYN_CORR_EXCL set");
		}

		// And clear the vector of observations
		vectorOfSourceObservations = null;
		vectorOfDestinationObservations = null;
	}
	
	/**
	 * <p>Computes the average Transfer Entropy for the previously supplied observations</p> 
	 * 
	 */
	public double computeAverageLocalOfObservations() throws Exception {
		
		double te = 0.0;
		if (debug) {
			MatrixUtils.printMatrix(System.out, destNextPastSourceVectors);
		}
		for (int b = 0; b < totalObservations; b++) {
			// Get the probability counts for the bth sample tuple
			int[] counts = getCounts(destNextPastSourceVectors[b], b);
			double logTerm = 0.0;
			double cont = 0.0;
			if (counts[NEXT_PAST_SOURCE] > 0) {
				logTerm = ((double) counts[NEXT_PAST_SOURCE] / (double) counts[PAST_SOURCE]) /
							((double) counts[NEXT_PAST] / (double) counts[PAST]);
				cont = Math.log(logTerm);
			}
			te += cont;
			if (debug) {
				System.out.println(b + ": " + destNextPastSourceVectors[b][0] + " (" + counts[NEXT_PAST_SOURCE] + " / " +
						counts[PAST_SOURCE] + ") / (" +
						counts[NEXT_PAST] + " / " + counts[PAST] + ") = " + 
						logTerm + " -> " + (cont/Math.log(2.0)) + " -> sum: " + (te/Math.log(2.0)));
			}
		}
		lastAverage = te / (double) totalObservations / Math.log(2.0);
		return lastAverage;
	}

	/**
	 * <p>Computes the average Transfer Entropy for the previously supplied observations,
	 *  using the Grassberger correction for the point count k: log_e(k) ~= digamma(k).</p>
	 * <p>Kaiser and Schreiber suggest (on p. 57) that for the TE
	 *   though the adverse correction of the bias correction is worse than the correction
	 *   itself (because the probabilities being multiplied/divided are not independent), 
	 *   so recommend not to use this method.
	 * </p>
	 * <p>It is implemented here for testing purposes only.</p>
	 * 
	 * @see Kaiser and Schreiber, Physica D 166 (2002) pp. 43-62
	 */
	public double computeAverageLocalOfObservationsWithCorrection() throws Exception {
		
		double te = 0.0;
		if (debug) {
			MatrixUtils.printMatrix(System.out, destNextPastSourceVectors);
		}
		for (int b = 0; b < totalObservations; b++) {
			int[] counts = getCounts(destNextPastSourceVectors[b], b);
			double cont = 0.0;
			if (counts[NEXT_PAST_SOURCE] > 0) {
				cont = MathsUtils.digamma(counts[NEXT_PAST_SOURCE]) -
						MathsUtils.digamma(counts[PAST_SOURCE]) -
						MathsUtils.digamma(counts[NEXT_PAST]) +
						MathsUtils.digamma(counts[PAST]);
			}
			te += cont;
			/*
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + (cont/Math.log(2.0)) + " -> sum: " + (te/Math.log(2.0)));
			}
			*/
		}
		// Average it, and convert results to bytes
		lastAverage = te / (double) totalObservations / Math.log(2.0);
		return lastAverage;
	}

	/**
	 * Computes the local transfer entropies for the previous supplied observations.
	 * 
	 * Where more than one time series has been added, the array
	 *  contains the local values for each tuple in the order in
	 *  which they were added.
	 * 
	 * If there was only a single time series added, the array
	 *  contains k zero values before the local values.
	 *  (This means the length of the return array is the same
	 *  as the length of the input time series).
	 * 
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		return computeLocalUsingPreviousObservations(null, null, true);
	}

	/**
	 * Comptues local transfer entropies for the given observations, using the previously supplied
	 * observations to compute the PDFs.
	 * I don't think it's such a good idea to do this for continuous variables (e.g. where
	 * one can get kernel estimates for probabilities of zero now) but I've implemented
	 * it anyway. I guess getting kernel estimates of zero here is no different than what
	 * can occur with dynamic correlation exclusion.
	 * 
	 * @param source
	 * @param destination
	 * @return
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(double[] source, double[] destination) throws Exception {
		return computeLocalUsingPreviousObservations(source, destination, false);
	}
	
	/**
	 * Returns the local TE at every time point.
	 * 
	 * @param source
	 * @param destination
	 * @param isPreviousObservations
	 * @return
	 * @throws Exception
	 */
	private double[] computeLocalUsingPreviousObservations(double[] source, double[] destination, boolean isPreviousObservations) throws Exception {
		
		double[][] newDestNextPastSourceVectors;

		if (isPreviousObservations) {
			// We've already computed the joint vectors for these observations
			newDestNextPastSourceVectors = destNextPastSourceVectors;
		} else {
			// We need to compute a new set of joint vectors
			newDestNextPastSourceVectors = makeJointVectorForNextPastSource(destination, source);
		}

		double te = 0.0;
		int numLocalObservations = newDestNextPastSourceVectors.length;
		double[] localTE;
		int offset = 0;
		if (isPreviousObservations && addedMoreThanOneObservationSet) {
			// We're returning the local values for a set of disjoint
			//  observations. So we don't add k zeros to the start
			localTE = new double[numLocalObservations];
			offset = 0;
		} else {
			localTE = new double[numLocalObservations + k];
			offset = k;
		}
		for (int b = 0; b < numLocalObservations; b++) {
			int[] counts = getCounts(newDestNextPastSourceVectors[b], isPreviousObservations ? b : -1);
			double logTerm = 0.0;
			double local = 0.0;
			if (counts[NEXT_PAST_SOURCE] > 0) {
				logTerm = ((double) counts[NEXT_PAST_SOURCE] / (double) counts[PAST_SOURCE]) /
							((double) counts[NEXT_PAST] / (double) counts[PAST]);
				local = Math.log(logTerm);
			}
			localTE[offset + b] = local;
			te += local;
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + (local/Math.log(2.0)) + " -> sum: " + (te/Math.log(2.0)));
			}
		}
		lastAverage = te / (double) numLocalObservations / Math.log(2.0);
		return localTE;
	}

	/**
	 * Compute the next, past and source values into a joint vector
	 * 
	 * @param destination
	 * @param source
	 * @return
	 */
	private double[][] makeJointVectorForNextPastSource(double[] destination, double[] source) {
		double[][] destNextPastSourceVectors = new double[destination.length - k][k + 2];
		for (int t = k; t < destination.length; t++) {
			for (int i = 0; i < k + 1; i++) {
				destNextPastSourceVectors[t - k][i] = destination[t - i];
			}
			destNextPastSourceVectors[t - k][k + 1] = source[t - 1];
		}
		return destNextPastSourceVectors;
	}
	
	/**
	 * Return the set of probability counts for the observation at time step 
	 * observationNumber checked against those stored in
	 * destNextPastSourceVectors.
	 * 
	 * @param observation
	 * @param observationNumber if &lt; 0 we don't try dynamic correlation exclusion
	 * @return
	 */
	protected int[] getCounts(double[] observation, int observationIndex) {
		int[] counts = new int[4];
		for (int b = 0; b < totalObservations; b++) {
			// If this is one of our samples, and we're doing dynamic correlation exclusion
			//  and it's within the window, don't count it
			if ((observationIndex >= 0) && dynCorrExcl &&  
					(Math.abs(b - observationIndex) < dynCorrExclTime)) {
				// This sample is within the dyn corr excl window,
				//  so don't count it
				continue;
			}

			if (!compareStateVectors(observation, b, 1, k)) {
				// Past vectors are different, no point continuing the comparison
				continue;
			}
			counts[PAST]++;
			boolean nextMatches = Math.abs(observation[0] - 
					destNextPastSourceVectors[b][0]) <= epsilons[0];
			boolean sourceMatches = Math.abs(observation[k + 1] - 
					destNextPastSourceVectors[b][k + 1]) <= epsilons[k + 1];
			counts[NEXT_PAST] += nextMatches ? 1 : 0;
			counts[PAST_SOURCE] += sourceMatches ? 1 : 0;
			counts[NEXT_PAST_SOURCE] += (nextMatches && sourceMatches) ? 1 : 0;
		}
		return counts;
	}
	
	/**
	 * Return whether the vectors at destNextPastSourceVectors[sampleIndex] and
	 * the observation are the same for length values from
	 * fromOffset
	 * 
	 * @param observation
	 * @param sampleIndex
	 * @param fromOffset
	 * @param length
	 * @return
	 */
	protected boolean compareStateVectors(double[] observation,
			int sampleIndex, int fromOffset, int length) {
		int i;
		for (i = 0; i < length; i++) {
			if (Math.abs(observation[fromOffset + i] - 
					destNextPastSourceVectors[sampleIndex][fromOffset + i]) > epsilons[fromOffset + i]) {
				// This pair can't be within the kernel
				return false;
			}
		}
		return true;
	}

	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		throw new RuntimeException("Not implemented in this calculator");
	}
	
	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception {
		throw new RuntimeException("Not implemented in this calculator");
	}

	/**
	 * <p>Computes the probability counts for each previously supplied observations</p>
	 * <p>Implemented primarily for debug purposes</p>
	 * 
	 */
	public KernelCount[][] computeMatchesForEachObservations(boolean giveListOfCorrelatedPoints) throws Exception {
		KernelCount[][] counts = new KernelCount[totalObservations][4];
		for (int b = 0; b < totalObservations; b++) {
			int[] intCounts = getCounts(destNextPastSourceVectors[b], b);
			counts[b] = new KernelCount[4];
			for (int i = 0; i < 4; i++) {
				int totalObsCount = 0;
				if (dynCorrExcl) {
					// Need to remove any observations that were *closer* than timeProximityForDynamicCorrelationExclusion
					int closeTimePointsToCompare = (b >= dynCorrExclTime) ?
							dynCorrExclTime - 1: b;
					closeTimePointsToCompare += (totalObservations - b >= dynCorrExclTime) ?
							dynCorrExclTime - 1: totalObservations - b - 1;
					closeTimePointsToCompare++; // Add one for comparison to self
					totalObsCount = totalObservations - closeTimePointsToCompare;
				} else {
					totalObsCount = totalObservations;
				}
				counts[b][i] = new KernelCount(intCounts[i], totalObsCount);
				// And ignore giveListOfCorrelatedPoints for now.
			}
		}
		return counts;
	}
}
