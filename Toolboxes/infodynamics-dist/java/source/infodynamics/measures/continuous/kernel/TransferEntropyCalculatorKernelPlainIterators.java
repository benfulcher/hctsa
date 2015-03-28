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
 * This class wraps the original observations in an iterator, to save on memory cost
 * of constructing joint vectors, but this results in a significant time cost
 * over {@link TransferEntropyCalculatorKernelPlain}.
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

/**
 * 
 * @author Joseph Lizier, joseph.lizier at gmail.com
 * @see For transfer entropy: Schreiber, PRL 85 (2) pp.461-464, 2000; http://dx.doi.org/10.1103/PhysRevLett.85.461
 * @see For local transfer entropy: Lizier et al, PRE 77, 026110, 2008; http://dx.doi.org/10.1103/PhysRevE.77.026110
 *
 */
public class TransferEntropyCalculatorKernelPlainIterators
	extends TransferEntropyCommon implements TransferEntropyCalculator {

	// Keep a record of which time series and time index in it
	//  each observation number corresponds to
	protected int[] timeSeriesIndex;
	protected int[] timeStepIndex;
	
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
	private double epsDest;
	private double epsSource;
	public static final String EPSILON_PROP_NAME = "EPSILON";
	
	/**
	 * Creates a new instance of the kernel-estimate style transfer entropy calculator
	 *
	 */
	public TransferEntropyCalculatorKernelPlainIterators() {
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
		super.initialise(k); // calls initialise();
	}

	/**
	 * Initialise using default or existing values for k and epsilon
	 */
	public void initialise() {
		timeSeriesIndex = null;
		timeStepIndex = null;
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
		super.setProperty(propertyName, propertyValue);
		if (propertyName.equalsIgnoreCase(EPSILON_PROP_NAME)) {
			epsilon = Double.parseDouble(propertyValue);
		} else if (propertyName.equalsIgnoreCase(NORMALISE_PROP_NAME)) {
			normalise = Boolean.parseBoolean(propertyValue);
		} else if (propertyName.equalsIgnoreCase(DYN_CORR_EXCL_TIME_NAME)) {
			dynCorrExclTime = Integer.parseInt(propertyValue);
			dynCorrExcl = (dynCorrExclTime > 0);
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

		// Now track which time series and step in it each observation comes from
		timeSeriesIndex = new int[totalObservations];
		timeStepIndex = new int[totalObservations];
		int obs = 0;
		int currentTimeSeriesIndex = 0;
		double sumDest = 0.0, sumSqDest = 0.0, sumSource = 0.0, sumSqSource = 0.0;
		Iterator<double[]> iterator = vectorOfSourceObservations.iterator();
		for (double[] destination : vectorOfDestinationObservations) {
			for (int t = k; t < destination.length; t++) {
				timeSeriesIndex[obs] = currentTimeSeriesIndex;
				timeStepIndex[obs++] = t;
			}
			currentTimeSeriesIndex++;
			if (normalise) {
				double[] source = iterator.next();
				// Track the sum and sum of squares of the dest and source
				// First add in the first k - 1 steps of history (dest only)
				for (int t = 0; t < k - 1; t++) {
					sumDest += destination[t];
					sumSqDest += destination[t] * destination[t];
				}
				// Now do all of the previous steps
				for (int t = k - 1; t < destination.length - 1; t++) {
					sumDest += destination[t];
					sumSqDest += destination[t] * destination[t];
					sumSource += source[t];
					sumSqSource += source[t] * source[t];
				}
				// Finally add the last destination step (dest only)
				sumDest += destination[destination.length - 1];
				sumSqDest += destination[destination.length - 1] * destination[destination.length - 1];
			}
		}
		
		if (normalise) {
			double meanDest = sumDest / (double) totalObservations;
			double meanSource = sumSource / (double) totalObservations;
			double meanSqDest = sumSqDest / (double) totalObservations;
			double meanSqSource = sumSqSource / (double) totalObservations;
			double stdDest = Math.sqrt(meanSqDest - meanDest*meanDest);
			double stdSource = Math.sqrt(meanSqSource - meanSource*meanSource);
			// Adjust epsilon for the dest and source
			epsDest = epsilon * stdDest;
			epsSource = epsilon * stdSource;
		} else {
			epsDest = epsilon;
			epsSource = epsilon;
		}
		
		// Store whether there was more than one observation set:
		addedMoreThanOneObservationSet = vectorOfDestinationObservations.size() > 1;
		if (addedMoreThanOneObservationSet && dynCorrExcl) {
			// This is fine for this calculator - see getCounts
		}
	}
	
	/**
	 * <p>Computes the average Transfer Entropy for the previously supplied observations</p> 
	 * 
	 */
	public double computeAverageLocalOfObservations() throws Exception {
		
		double te = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			int timeSeries = timeSeriesIndex[b];
			double[] source = vectorOfSourceObservations.elementAt(timeSeries);
			double[] dest = vectorOfDestinationObservations.elementAt(timeSeries);
			int[] counts = getCounts(source, dest, timeStepIndex[b], timeSeries);
			double logTerm = 0.0;
			double cont = 0.0;
			if (counts[NEXT_PAST_SOURCE] > 0) {
				logTerm = ((double) counts[NEXT_PAST_SOURCE] / (double) counts[PAST_SOURCE]) /
							((double) counts[NEXT_PAST] / (double) counts[PAST]);
				cont = Math.log(logTerm);
			}
			te += cont;
			if (debug) {
				System.out.println(b + ": " + dest[timeStepIndex[b]] + " (" + counts[NEXT_PAST_SOURCE] + " / " +
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
	 * <p>Kaiser and Schreiber, Physica D 166 (2002) pp. 43-62 suggest (on p. 57) that for the TE
	 *   though the adverse correction of the bias correction is worse than the correction
	 *   itself (because the probabilities being multiplied/divided are not independent), 
	 *   so recommend not to use this method.
	 * </p>
	 * <p>It is implemented here for testing purposes only.</p>
	 * 
	 */
	public double computeAverageLocalOfObservationsWithCorrection() throws Exception {
		
		double te = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			int timeSeries = timeSeriesIndex[b];
			double[] source = vectorOfSourceObservations.elementAt(timeSeries);
			double[] dest = vectorOfDestinationObservations.elementAt(timeSeries);
			int[] counts = getCounts(source, dest, timeStepIndex[b], timeSeries);
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
		
		int numLocalObservations;
		if (isPreviousObservations) {
			numLocalObservations = totalObservations;
		} else {
			numLocalObservations = destination.length - k;
		}

		double te = 0.0;
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
			int timeSeries = -1;
			int timeStep = b;
			if (isPreviousObservations) {
				timeSeries = timeSeriesIndex[b];
				timeStep = timeStepIndex[b];
				source = vectorOfSourceObservations.elementAt(timeSeries);
				destination = vectorOfDestinationObservations.elementAt(timeSeries);
			}
			int[] counts = getCounts(source, destination, timeStep, timeSeries);
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
	 * Return the set of probability counts for the observation at time step 
	 * timeStep of the given source and destination checked against those stored in
	 * the supplied samples.
	 * timeSeries tells us which sample these were from, in case we want to do
	 * dynamic correlation exclusion. If it is -1, then they are not from one of the samples.
	 * 
	 * @param observation
	 * @param observationNumber
	 * @return
	 */
	protected int[] getCounts(double[] source, double[] dest, int timeStep, int timeSeries) {
		int[] counts = new int[4];
		for (int b = 0; b < totalObservations; b++) {
			// Pull out the current source and dest samples to compare to
			double[] sampleSource = vectorOfSourceObservations.elementAt(timeSeriesIndex[b]);
			double[] sampleDest = vectorOfDestinationObservations.elementAt(timeSeriesIndex[b]);
			int t = timeStepIndex[b];
			// If this is one of our samples, and we're doing dynamic correlation exclusion,
			//  and they're from the same sample time series,
			//  and it's within the window ==> don't count it
			if ((timeSeries >= 0) && dynCorrExcl &&
					(timeSeriesIndex[b] == timeSeries) &&
					(Math.abs(t - timeStep) < dynCorrExclTime)) {
				// This sample is within the dyn corr excl window,
				//  so don't count it
				continue;
			}

			if (!compareHistory(dest, timeStep, sampleDest, t)) {
				// Past vectors are different, no point continuing the comparison
				continue;
			}
			counts[PAST]++;
			boolean nextMatches = Math.abs(dest[timeStep] - sampleDest[t]) <= epsDest;
			boolean sourceMatches = Math.abs(source[timeStep - 1] - sampleSource[t - 1])
										<= epsSource;
			counts[NEXT_PAST] += nextMatches ? 1 : 0;
			counts[PAST_SOURCE] += sourceMatches ? 1 : 0;
			counts[NEXT_PAST_SOURCE] += (nextMatches && sourceMatches) ? 1 : 0;
		}
		return counts;
	}
	
	/**
	 * Return whether the past k time steps (prior to t1 and t2) of the
	 *  two time series series1 and series2 match within epsDest at each point
	 * 
	 * @param series1
	 * @param t1
	 * @param series2
	 * @param t2
	 * @return
	 */
	protected boolean compareHistory(double[] series1, int t1, double[] series2, int t2) {
		for (int i = 1; i <= k; i++) {
			if (Math.abs(series1[t1 - i] - series2[t2 - i]) > epsDest) {
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
			int timeSeries = timeSeriesIndex[b];
			double[] source = vectorOfSourceObservations.elementAt(timeSeries);
			double[] dest = vectorOfDestinationObservations.elementAt(timeSeries);
			int[] intCounts = getCounts(source, dest, timeStepIndex[b], timeSeries);
			counts[b] = new KernelCount[4];
			for (int i = 0; i < 4; i++) {
				int totalObsCount = 0;
				if (dynCorrExcl) {
					// Need to remove any observations that were in the same
					//  sample time series and *closer* than dynCorrExclTime
					int observationsInSameTimeSeries = dest.length;
					int obsIndexInItsTimeSeries = timeStepIndex[b] - k;
					int closeTimePointsToCompare = (obsIndexInItsTimeSeries >= dynCorrExclTime) ?
							dynCorrExclTime - 1: obsIndexInItsTimeSeries;
					closeTimePointsToCompare += (observationsInSameTimeSeries - obsIndexInItsTimeSeries >= dynCorrExclTime) ?
							dynCorrExclTime - 1: observationsInSameTimeSeries - obsIndexInItsTimeSeries - 1;
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
