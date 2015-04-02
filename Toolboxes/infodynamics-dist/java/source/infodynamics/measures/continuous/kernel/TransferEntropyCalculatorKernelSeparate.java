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
import infodynamics.measures.continuous.kernel.KernelEstimatorMultiVariate;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

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
 * <p>***** This implementation is simple in terms of using separate
 * kernel estimators for each probability calculation required.
 * The main TE box-kernel class, {@link TransferEntropyCalculatorKernel} is significantly
 * more efficient.
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
 * <p>
 * TODO Implement dynamic correlation exclusion with multiple observation sets. (see the
 *  way this is done in Plain calculator).
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
public class TransferEntropyCalculatorKernelSeparate
	extends TransferEntropyCommon implements TransferEntropyCalculator {

	protected KernelEstimatorMultiVariate mvkeDestinationPast = null;
	protected KernelEstimatorMultiVariate mvkeDestinationNextPast = null;
	protected KernelEstimatorMultiVariate mvkeDestinationPastSource = null;
	protected KernelEstimatorMultiVariate mvkeDestinationNextPastSource = null;

	// Keep joint vectors so we don't need to regenerate them
	protected double[][] destPastVectors;
	protected double[][] destNextPastVectors;
	protected double[][] destPastSourceVectors;
	protected double[][] destNextPastSourceVectors;
	
	private boolean normalise = true;
	public static final String NORMALISE_PROP_NAME = "NORMALISE";
	
	private boolean dynCorrExcl = false;
	private int dynCorrExclTime = 100;
	public static final String DYN_CORR_EXCL_TIME_NAME = "DYN_CORR_EXCL";
	
	private boolean forceCompareToAll = false;
	public static final String FORCE_KERNEL_COMPARE_TO_ALL = "FORCE_KERNEL_COMPARE_TO_ALL";
	
	/**
	 * Default value for epsilon
	 */
	public static final double DEFAULT_EPSILON = 0.25;
	/**
	 * Kernel width
	 */
	private double epsilon = DEFAULT_EPSILON;
	public static final String EPSILON_PROP_NAME = "EPSILON";
	
	/**
	 * Creates a new instance of the kernel-estimate style transfer entropy calculator
	 *
	 */
	public TransferEntropyCalculatorKernelSeparate() {
		super();
		mvkeDestinationPast = new KernelEstimatorMultiVariate();
		mvkeDestinationNextPast = new KernelEstimatorMultiVariate();
		mvkeDestinationPastSource = new KernelEstimatorMultiVariate();
		mvkeDestinationNextPastSource = new KernelEstimatorMultiVariate();
		mvkeDestinationPast.setNormalise(normalise);
		mvkeDestinationNextPast.setNormalise(normalise);
		mvkeDestinationPastSource.setNormalise(normalise);
		mvkeDestinationNextPastSource.setNormalise(normalise);
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
	public void initialise() throws Exception {
		mvkeDestinationPast.initialise(k, epsilon);
		mvkeDestinationNextPast.initialise(k+1, epsilon);
		mvkeDestinationPastSource.initialise(k+1, epsilon);
		mvkeDestinationNextPastSource.initialise(k+2, epsilon);
		destPastVectors = null;
		destNextPastVectors = null;
		destPastSourceVectors = null;
		destNextPastSourceVectors = null;
	}
	
	/**
	 * Set properties for the transfer entropy calculator.
	 * These can include:
	 * <ul>
	 * 		<li>K_PROP_NAME</li>
	 * 		<li>EPSILON_PROP_NAME</li>
	 * 		<li>NORMALISE_PROP_NAME</li>
	 * 		<li>DYN_CORR_EXCL_TIME_NAME</li>
	 * 		<li>FORCE_KERNEL_COMPARE_TO_ALL</li>
	 * </ul> 
	 * 
	 * @param propertyName
	 * @param propertyValue
	 * @throws Exception 
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		super.setProperty(propertyName, propertyValue);
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(EPSILON_PROP_NAME)) {
			epsilon = Double.parseDouble(propertyValue);
		} else if (propertyName.equalsIgnoreCase(NORMALISE_PROP_NAME)) {
			normalise = Boolean.parseBoolean(propertyValue);
			mvkeDestinationPast.setNormalise(normalise);
			mvkeDestinationNextPast.setNormalise(normalise);
			mvkeDestinationPastSource.setNormalise(normalise);
			mvkeDestinationNextPastSource.setNormalise(normalise);
		} else if (propertyName.equalsIgnoreCase(DYN_CORR_EXCL_TIME_NAME)) {
			dynCorrExclTime = Integer.parseInt(propertyValue);
			dynCorrExcl = (dynCorrExclTime > 0);
			if (dynCorrExcl) {
				mvkeDestinationPast.setDynamicCorrelationExclusion(dynCorrExclTime);
				mvkeDestinationNextPast.setDynamicCorrelationExclusion(dynCorrExclTime);
				mvkeDestinationPastSource.setDynamicCorrelationExclusion(dynCorrExclTime);
				mvkeDestinationNextPastSource.setDynamicCorrelationExclusion(dynCorrExclTime);
			} else {
				mvkeDestinationPast.clearDynamicCorrelationExclusion();
				mvkeDestinationNextPast.clearDynamicCorrelationExclusion();
				mvkeDestinationPastSource.clearDynamicCorrelationExclusion();
				mvkeDestinationNextPastSource.clearDynamicCorrelationExclusion();
			}
		} else if (propertyName.equalsIgnoreCase(FORCE_KERNEL_COMPARE_TO_ALL)) {
			forceCompareToAll = Boolean.parseBoolean(propertyValue);
			mvkeDestinationPast.setForceCompareToAll(forceCompareToAll);
			mvkeDestinationNextPast.setForceCompareToAll(forceCompareToAll);
			mvkeDestinationPastSource.setForceCompareToAll(forceCompareToAll);
			mvkeDestinationNextPastSource.setForceCompareToAll(forceCompareToAll);
		} else {
			// No property was set
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println("Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	/*
	 * Old implementation of set observations - now we defer to the super class
	 *  and let it call startAdd, Add and FinaliseAdd.
	 * 
	 * ------------------------
	 * 
	 * Sets the observations to compute the PDFs from.
	 * Cannot be called in conjunction with start/add/finaliseAddObservations
	 * 
	 * @param source observations for the source variable
	 * @param destination observations for the destination variable
	 * 
	public void setObservations(double[] source, double[] destination) {
		totalObservations = destination.length - k;
		// We're just going to be wasteful with memory for the moment since it's easier to
		//  implement for now:
		// Construct joint vectors past
		destPastVectors = makeJointVectorForPast(destination);
		mvkeDestinationPast.setObservations(destPastVectors);
		// Construct joint vectors for next and past
		destNextPastVectors = makeJointVectorForNextPast(destination);
		mvkeDestinationNextPast.setObservations(destNextPastVectors);
		// Construct joint vectors for past and source
		destPastSourceVectors = makeJointVectorForPastSource(destination, source);
		mvkeDestinationPastSource.setObservations(destPastSourceVectors);
		// Construct joint vectors for next and past and source
		destNextPastSourceVectors = makeJointVectorForNextPastSource(destination, source);
		mvkeDestinationNextPastSource.setObservations(destNextPastSourceVectors);
	}
	*/

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
		destPastVectors = new double[totalObservations][k];
		destNextPastVectors = new double[totalObservations][k + 1];
		destPastSourceVectors = new double[totalObservations][k + 1];
		destNextPastSourceVectors = new double[totalObservations][k + 2];
		
		// Construct the joint vectors from the given observations
		int startObservation = 0;
		Iterator<double[]> iterator = vectorOfDestinationObservations.iterator();
		for (double[] source : vectorOfSourceObservations) {
			double[] destination = iterator.next();
			double[][] currentDestPastVectors = makeJointVectorForPast(destination);
			MatrixUtils.arrayCopy(currentDestPastVectors, 0, 0,
					destPastVectors, startObservation, 0, currentDestPastVectors.length, k);
			double[][] currentDestNextPastVectors = makeJointVectorForNextPast(destination);
			MatrixUtils.arrayCopy(currentDestNextPastVectors, 0, 0,
					destNextPastVectors, startObservation, 0, currentDestNextPastVectors.length, k + 1);
			double[][] currentDestPastSourceVectors = makeJointVectorForPastSource(destination, source);
			MatrixUtils.arrayCopy(currentDestPastSourceVectors, 0, 0,
					destPastSourceVectors, startObservation, 0, currentDestPastSourceVectors.length, k + 1);
			double[][] currentDestNextPastSourceVectors = makeJointVectorForNextPastSource(destination, source);
			MatrixUtils.arrayCopy(currentDestNextPastSourceVectors, 0, 0,
					destNextPastSourceVectors, startObservation, 0, currentDestNextPastSourceVectors.length, k + 2);
			startObservation += destination.length - k;
		}
		
		// Now set the joint vectors in the kernel estimators
		mvkeDestinationPast.setObservations(destPastVectors);
		mvkeDestinationNextPast.setObservations(destNextPastVectors);
		mvkeDestinationPastSource.setObservations(destPastSourceVectors);
		mvkeDestinationNextPastSource.setObservations(destNextPastSourceVectors);

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
			double countPast = mvkeDestinationPast.getCount(destPastVectors[b], b);
			double countNextPast = mvkeDestinationNextPast.getCount(destNextPastVectors[b], b);
			double countPastSource = mvkeDestinationPastSource.getCount(destPastSourceVectors[b], b);
			double countNextPastSource = mvkeDestinationNextPastSource.getCount(destNextPastSourceVectors[b], b);
			double logTerm = 0.0;
			double cont = 0.0;
			if (countNextPastSource > 0) {
				logTerm = (countNextPastSource / countPastSource) / (countNextPast / countPast);
				cont = Math.log(logTerm);
			}
			te += cont;
			if (debug) {
				System.out.println(b + ": " + destPastVectors[b][0] + " (" + countNextPastSource + " / " + countPastSource + ") / (" +
						countNextPast + " / " + countPast + ") = " + 
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
		if (debug) {
			MatrixUtils.printMatrix(System.out, destNextPastSourceVectors);
		}
		for (int b = 0; b < totalObservations; b++) {
			int countPast = mvkeDestinationPast.getCount(destPastVectors[b], b);
			int countNextPast = mvkeDestinationNextPast.getCount(destNextPastVectors[b], b);
			int countPastSource = mvkeDestinationPastSource.getCount(destPastSourceVectors[b], b);
			int countNextPastSource = mvkeDestinationNextPastSource.getCount(destNextPastSourceVectors[b], b);
			double cont = 0.0;
			if (countNextPastSource > 0) {
				cont = MathsUtils.digamma(countNextPastSource) - 
						MathsUtils.digamma(countPastSource) -
						MathsUtils.digamma(countNextPast) +
						MathsUtils.digamma(countPast);
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

	/*
	 * TEST ONLY - THIS SHOULD not be used:
	 * 
	 * This is the way I previsouly computed TE in octave code.
	 * This method of estimating p(i_n, i_n+1, j_n) should not be correct, because
	 *  there is no need to multiply by the prob of observation.
	 * I'm just writing it in here to test .. 
	 * 
	 * @return
	 * @throws Exception
	private double computeAverageLocalOfObservationsWithMultiplier() throws Exception {

		double te = 0.0;
		double sumJointProb = 0;
		for (int b = 0; b < totalObservations; b++) {
			double probPast = mvkeDestinationPast.getProbability(destPastVectors[b], b);
			double probNextPast = mvkeDestinationNextPast.getProbability(destNextPastVectors[b], b);
			double probPastSource = mvkeDestinationPastSource.getProbability(destPastSourceVectors[b], b);
			double probNextPastSource = mvkeDestinationNextPastSource.getProbability(destNextPastSourceVectors[b], b);
			double logTerm = 0.0;
			double cont = 0.0;
			if (probNextPastSource > ZERO_COMPARATOR) {
				logTerm = (probNextPastSource / probPastSource) / (probNextPast / probPast);
				cont = Math.log(logTerm);
			}
			sumJointProb += probNextPastSource;
			te += cont * probNextPastSource;
		}
		lastAverage = te / sumJointProb / Math.log(2.0);
		return lastAverage;
	}
	*/
	
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
		
		double[][] newDestPastVectors;
		double[][] newDestNextPastVectors;
		double[][] newDestPastSourceVectors;
		double[][] newDestNextPastSourceVectors;

		if (isPreviousObservations) {
			// We've already computed the joint vectors for these observations
			newDestPastVectors = destPastVectors;
			newDestNextPastVectors = destNextPastVectors;
			newDestPastSourceVectors = destPastSourceVectors;
			newDestNextPastSourceVectors = destNextPastSourceVectors;
		} else {
			// We need to compute a new set of joint vectors
			newDestPastVectors = makeJointVectorForPast(destination);
			newDestNextPastVectors = makeJointVectorForNextPast(destination);
			newDestPastSourceVectors = makeJointVectorForPastSource(destination, source);
			newDestNextPastSourceVectors = makeJointVectorForNextPastSource(destination, source);
		}

		double te = 0.0;
		int numLocalObservations = newDestPastVectors.length;
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
		double countPast, countNextPast, countPastSource, countNextPastSource;
		for (int b = 0; b < numLocalObservations; b++) {
			if (isPreviousObservations) {
				countPast = mvkeDestinationPast.getCount(newDestPastVectors[b], b);
				countNextPast = mvkeDestinationNextPast.getCount(newDestNextPastVectors[b], b);
				countPastSource = mvkeDestinationPastSource.getCount(newDestPastSourceVectors[b], b);
				countNextPastSource = mvkeDestinationNextPastSource.getCount(newDestNextPastSourceVectors[b], b);
			} else {
				countPast = mvkeDestinationPast.getCount(newDestPastVectors[b]);
				countNextPast = mvkeDestinationNextPast.getCount(newDestNextPastVectors[b]);
				countPastSource = mvkeDestinationPastSource.getCount(newDestPastSourceVectors[b]);
				countNextPastSource = mvkeDestinationNextPastSource.getCount(newDestNextPastSourceVectors[b]);
			}
			double logTerm = 0.0;
			double local = 0.0;
			if (countNextPastSource > 0) {
				logTerm = (countNextPastSource / countPastSource) / (countNextPast / countPast);
				local = Math.log(logTerm) / Math.log(2.0);
			}
			localTE[offset + b] = local;
			te += local;
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + local + " -> sum: " + te);
			}
		}
		lastAverage = te / (double) numLocalObservations;
		return localTE;
	}

	/**
	 * Combine the past and source values into a joint vector
	 * 
	 * @param destination
	 * @param source
	 * @return
	 */
	private double[][] makeJointVectorForPastSource(double[] destination, double[] source) {
		double[][] destPastSourceVectors = new double[destination.length - k][k + 1];
		for (int t = k; t < destination.length; t++) {
			for (int i = 0; i < k; i++) {
				destPastSourceVectors[t - k][i] = destination[t - i - 1];
			}
			destPastSourceVectors[t - k][k] = source[t - 1];
		}
		return destPastSourceVectors;
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
	 * Compute the significance of obtaining the given average TE from the given observations
	 * 
	 * 	This is as per Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128.
	 *
	 * Basically, we shuffle the source observations against the destination tuples.
	 * This keeps the marginal PDFs the same (including the entropy rate of the destination)
	 *  but destroys any correlation between the source and state change of the destination.
	 * 
	 * @param numPermutationsToCheck number of new orderings of the source values to compare against
	 * @return
	 */
	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(totalObservations, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}
	
	/**
	 * As per {@link computeSignificance(int) computeSignificance()} but supplies
	 *  the re-orderings of the observations of the source variables.
	 *  
	 * 
	 * @param newOrderings first index is permutation number, i.e. newOrderings[i]
	 * 		is an array of 1 permutation of 0..n-1, where there were n observations.
	 * @return
	 * @throws Exception
	 */
	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception {
		
		int numPermutationsToCheck = newOrderings.length;
				
		double actualTE = computeAverageLocalOfObservations();
		
		// Save the relevant source observations here:
		double[] originalSourceValuesInJoint = MatrixUtils.selectColumn(destPastSourceVectors, k);
		
		int countWhereTeIsMoreSignificantThanOriginal = 0;
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		for (int p = 0; p < numPermutationsToCheck; p++) {
			// Check that the length of the reorderings is OK
			if (newOrderings[p].length != totalObservations) {
				throw new Exception("Length " + newOrderings[p].length +
						" of reordering " + p + " in newOrderings does not " +
						" match the number of observations " + totalObservations);
			}

			// Generate a new re-ordered data set for the source in the destPastSourceVectors 
			//  and destNextPastSourceVectors vectors
			MatrixUtils.reorderVectorIntoMatrix(originalSourceValuesInJoint, newOrderings[p],
					destPastSourceVectors, k);
			MatrixUtils.reorderVectorIntoMatrix(originalSourceValuesInJoint, newOrderings[p],
					destNextPastSourceVectors, k+1);
			
			// Make the equivalent operations of intialise
			mvkeDestinationPastSource.initialise(k+1, epsilon);
			mvkeDestinationNextPastSource.initialise(k+2, epsilon);
			// Make the equivalent operations of setObservations:
			mvkeDestinationPastSource.setObservations(destPastSourceVectors);
			mvkeDestinationNextPastSource.setObservations(destNextPastSourceVectors);
			// And get a TE value for this realisation:
			double newTe = computeAverageLocalOfObservations();
			measDistribution.distribution[p] = newTe;
			if (newTe >= actualTE) {
				countWhereTeIsMoreSignificantThanOriginal++;
			}
		}
		
		// Restore the local variables:
		lastAverage = actualTE;
		// Restore the source observations in the joint vectors
		MatrixUtils.copyIntoColumn(destPastSourceVectors, k, originalSourceValuesInJoint);
		MatrixUtils.copyIntoColumn(destNextPastSourceVectors, k+1, originalSourceValuesInJoint);
		// And set the mulit-variate kernel estimators back to their previous state
		mvkeDestinationPastSource.initialise(k+1, epsilon);
		mvkeDestinationNextPastSource.initialise(k+2, epsilon);
		mvkeDestinationPastSource.setObservations(destPastSourceVectors);
		mvkeDestinationNextPastSource.setObservations(destNextPastSourceVectors);
		
		// And return the significance
		measDistribution.pValue = (double) countWhereTeIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualTE;
		return measDistribution;
	}

	/**
	 * <p>Computes the probability counts for each previously supplied observations</p>
	 * <p>Implemented primarily for debug purposes</p>
	 * 
	 */
	public KernelCount[][] computeMatchesForEachObservations(boolean giveListOfCorrelatedPoints) throws Exception {
		
		KernelCount[][] counts = new KernelCount[totalObservations][4];
		for (int b = 0; b < totalObservations; b++) {
			counts[b][0] = mvkeDestinationPast.
				getCompleteKernelCount(destPastVectors[b], b, giveListOfCorrelatedPoints);
			counts[b][1] = mvkeDestinationNextPast.
				getCompleteKernelCount(destNextPastVectors[b], b, giveListOfCorrelatedPoints);
			counts[b][2] = mvkeDestinationPastSource.
				getCompleteKernelCount(destPastSourceVectors[b], b, giveListOfCorrelatedPoints);
			counts[b][3] = mvkeDestinationNextPastSource.
				getCompleteKernelCount(destNextPastSourceVectors[b], b, giveListOfCorrelatedPoints);
		}
		return counts;
	}

	public void setDebug(boolean debug) {
		super.setDebug(debug);
		mvkeDestinationPast.setDebug(debug);
		mvkeDestinationNextPast.setDebug(debug);
		mvkeDestinationPastSource.setDebug(debug);
		mvkeDestinationNextPastSource.setDebug(debug);
	}
}
