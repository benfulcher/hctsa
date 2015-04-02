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

import infodynamics.measures.continuous.TransferEntropyCalculatorMultiVariate;
import infodynamics.measures.continuous.TransferEntropyCommon;
import infodynamics.measures.continuous.kernel.TransferEntropyKernelCounts;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

import java.util.Iterator;
import java.util.Vector;

/**
 * <p>Computes the differential transfer entropy (TE) between two <b>multivariate</b>
 *  <code>double[][]</code> time-series of observations
 *  using box-kernel estimation.
 *  See Schreiber below for the definition of transfer entropy,
 *  Lizier et al. (2011) for the extension to multivariate source
 *  and destination, and
 *  and Lizier et al. (2008) for the definition of local transfer entropy.
 *  For details on box-kernel estimation, see Kantz and Schreiber (below).</p>
 *  
 * <p>TE was defined by Schreiber (below).
 *  This estimator is realised here by extending
 *  {@link TransferEntropyCommon}.</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link TransferEntropyCalculatorMultiVariate},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #TransferEntropyCalculatorMultiVariateKernel()}.</li>
 *  <li>Further properties are available, see {@link #setProperty(String, String)};</li>
 *  <li>Additional {@link #initialise(int)}, {@link #initialise(int, double)} options;</li>
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
 *  <li>J.T. Lizier, J. Heinzle, A. Horstmann, J.-D. Haynes, M. Prokopenko,
 *  <a href="http://dx.doi.org/10.1007/s10827-010-0271-2">
 *  "Multivariate information-theoretic measures reveal directed information
 *  structure and task relevant changes in fMRI connectivity"</a>,
 *  Journal of Computational Neuroscience, vol. 30, pp. 85-107, 2011.</li>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href="http://dx.doi.org/10.1103/PhysRevE.77.026110">
 *  "Local information transfer as a spatiotemporal filter for complex systems"</a>
 *  Physical Review E 77, 026110, 2008.</li>
 *  <li>H. Kantz and T. Schreiber, "Nonlinear Time Series Analysis"
 *  (Cambridge University Press, Cambridge, MA, 1997).</li>
 * </ul>
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class TransferEntropyCalculatorMultiVariateKernel
	extends TransferEntropyCommon implements TransferEntropyCalculatorMultiVariate {

	protected KernelEstimatorTransferEntropyMultiVariate teKernelEstimator = null;
	// Keep a kernel estimator for the next state, in case we wish to compute
	//  Active info storage also:
	protected KernelEstimatorMultiVariate nextStateKernelEstimator = null;

	/**
	 * Storage for source observations for addObservsations
	 */
	protected Vector<double[][]> vectorOfJointSourceObservations;
	/**
	 * Storage for destination observations for addObservsations
	 */
	protected Vector<double[][]> vectorOfJointDestinationObservations;

	// Keep joint vectors so we don't need to regenerate them
	protected double[][] destPastVectors;
	protected double[][] destNextVectors;
	protected double[][] sourceVectors;
	
	protected int destDimensions = 1;
	protected int sourceDimensions = 1;
	
	// Store the local conditional probability of next on past state as 
	//  computed during a local TE computation, in case the caller wants to
	//  compute the local active info storage next.
	protected double[] localProbNextCondPast;
	
	protected boolean normalise = true;
	public static final String NORMALISE_PROP_NAME = "NORMALISE";
	
	protected boolean dynCorrExcl = false;
	protected int dynCorrExclTime = 100;
	public static final String DYN_CORR_EXCL_TIME_NAME = "DYN_CORR_EXCL";
	
	protected boolean forceCompareToAll = false;
	public static final String FORCE_KERNEL_COMPARE_TO_ALL = "FORCE_KERNEL_COMPARE_TO_ALL";
	
	/**
	 * Default value for kernel width
	 */
	public static final double DEFAULT_KERNEL_WIDTH = 0.25;
	/**
	 * Kernel width
	 */
	protected double kernelWidth = DEFAULT_KERNEL_WIDTH;
	/**
	 * Property name for the kernel width
	 */
	public static final String KERNEL_WIDTH_PROP_NAME = "KERNEL_WIDTH";
	/**
	 * Legacy property name for the kernel width
	 */
	public static final String EPSILON_PROP_NAME = "EPSILON";
	
	/**
	 * Creates a new instance of the kernel-estimate style transfer entropy calculator
	 *
	 */
	public TransferEntropyCalculatorMultiVariateKernel() {
		super();
		teKernelEstimator = new KernelEstimatorTransferEntropyMultiVariate();
		teKernelEstimator.setNormalise(normalise);
		nextStateKernelEstimator = new KernelEstimatorMultiVariate();
		nextStateKernelEstimator.setNormalise(normalise);
	}

	@Override
	public void initialise(int k) throws Exception {
		initialise(k, kernelWidth);
	}
	
	/**
	 * Initialise the calculator for (re-)use, with a specific 
	 * embedded destination history length and kernel width,
	 * and existing (or default) values of other parameters,.
	 * Clears an PDFs of previously supplied observations.
	 *
	 * @param k destination embedded history length
	 * @param kernelWidth if {@link #NORMALISE_PROP_NAME} property has
	 *  been set, then this kernel width corresponds to the number of
	 *  standard deviations from the mean (otherwise it is an absolute value)
	 */
	public void initialise(int k, double kernelWidth) throws Exception {
		this.kernelWidth = kernelWidth;
		initialise(k, 1, 1); // assume 1 dimension in source and dest
	}

	@Override
	public void initialise(int k, int sourceDimensions, int destDimensions) throws Exception {
		this.destDimensions = destDimensions;
		this.sourceDimensions = sourceDimensions;
		super.initialise(k); // calls initialise();
	}
	
	@Override
	public void initialise(int sourceDimensions, int destDimensions) throws Exception {
		this.destDimensions = destDimensions;
		this.sourceDimensions = sourceDimensions;
		super.initialise(k); // calls initialise();
	}

	@Override
	public void initialise() {
		teKernelEstimator.initialise(k * destDimensions,
				sourceDimensions, kernelWidth, kernelWidth);
		nextStateKernelEstimator.initialise(destDimensions, kernelWidth);
		destPastVectors = null;
		destNextVectors = null;
		sourceVectors = null;
		localProbNextCondPast = null;
	}
	
	/**
	 * <p>Set properties for the kernel TE calculator.
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
	 * 		<li>{@link #NORMALISE_PROP_NAME} -- whether to normalise the incoming variables 
	 * 			to mean 0, standard deviation 1, or not (default false). Sets {@link #normalise}.</li>
	 * 		<li>{@link #DYN_CORR_EXCL_TIME_NAME} -- a dynamics exclusion time window (see Kantz and Schreiber),
	 * 			default is 0 which means no dynamic exclusion window.</li>
	 * 		<li>{@link #FORCE_KERNEL_COMPARE_TO_ALL} -- whether to force the underlying kernel estimators to compare
	 *  		each data point to each other (or else allow it to use optimisations).</li>
	 *  	<li>any valid properties for {@link TransferEntropyCommon#setProperty(String, String)}.</li>
	 * </ul>
	 * </p>
	 * 
	 * <p>Note that dynamic correlation exclusion (set with {@link #DYN_CORR_EXCL_TIME_NAME})
	 *  may have unexpected results if multiple
	 *  observation sets have been added. This is because multiple observation sets
	 *  are treated as though they are from a single time series, so observations from
	 *  near the end of observation set i will be excluded from comparison to 
	 *  observations near the beginning of observation set (i+1). 
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(KERNEL_WIDTH_PROP_NAME) ||
				propertyName.equalsIgnoreCase(EPSILON_PROP_NAME)) {
			kernelWidth = Double.parseDouble(propertyValue);
		} else if (propertyName.equalsIgnoreCase(NORMALISE_PROP_NAME)) {
			normalise = Boolean.parseBoolean(propertyValue);
			teKernelEstimator.setNormalise(normalise);
			nextStateKernelEstimator.setNormalise(normalise);
		} else if (propertyName.equalsIgnoreCase(DYN_CORR_EXCL_TIME_NAME)) {
			dynCorrExclTime = Integer.parseInt(propertyValue);
			dynCorrExcl = (dynCorrExclTime > 0);
			if (dynCorrExcl) {
				teKernelEstimator.setDynamicCorrelationExclusion(dynCorrExclTime);
				nextStateKernelEstimator.setDynamicCorrelationExclusion(dynCorrExclTime);
			} else {
				teKernelEstimator.clearDynamicCorrelationExclusion();
				nextStateKernelEstimator.clearDynamicCorrelationExclusion();
			}
		} else if (propertyName.equalsIgnoreCase(FORCE_KERNEL_COMPARE_TO_ALL)) {
			forceCompareToAll = Boolean.parseBoolean(propertyValue);
			teKernelEstimator.setForceCompareToAll(forceCompareToAll);
			nextStateKernelEstimator.setForceCompareToAll(forceCompareToAll);
		} else {
			// No property was set
			propertySet = false;
			// try the superclass:
			super.setProperty(propertyName, propertyValue);
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	@Override
	public void setObservations(double[][] source, double[][] destination) throws Exception {
		startAddObservations();
		addObservations(source, destination);
		finaliseAddObservations();
	}
	
	@Override
	public void setObservations(double[][] source, double[][] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception {
		Vector<int[]> startAndEndTimePairs = computeStartAndEndTimePairs(sourceValid, destValid);
		
		// We've found the set of start and end times for this pair
		startAddObservations();
		for (int[] timePair : startAndEndTimePairs) {
			int startTime = timePair[0];
			int endTime = timePair[1];
			addObservations(source, destination, startTime, endTime - startTime + 1);
		}
		finaliseAddObservations();
	}

	@Override
	public void setObservations(double[][] source, double[][] destination,
			boolean[][] sourceValid, boolean[][] destValid) throws Exception {
		
		boolean[] jointSourceValid = MatrixUtils.andRows(sourceValid);
		boolean[] jointDestValid = MatrixUtils.andRows(destValid);
		setObservations(source, destination, jointSourceValid, jointDestValid);
	}

	@Override
	public void startAddObservations() {
		vectorOfJointSourceObservations = new Vector<double[][]>();
		vectorOfJointDestinationObservations = new Vector<double[][]>();
	}

	/**
	 * Add observations of a single-dimensional source and destination pair.
	 * This call is only allowed if the source and destination
	 * dimensions are set to 1.
	 * 
	 * {@inheritDoc} 
	 */
	@Override
	public void addObservations(double[] source, double[] destination) throws Exception {
		double[][] sourceMatrix = new double[source.length][1];
		MatrixUtils.copyIntoColumn(sourceMatrix, 0, source);
		double[][] destMatrix = new double[destination.length][1];
		MatrixUtils.copyIntoColumn(destMatrix, 0, destination);
		addObservations(sourceMatrix, destMatrix);
	}

	/**
	 * Add sub-series of observations of a single-dimensional source and destination pair.
	 * This call is only allowed if the source and destination
	 * dimensions are set to 1.
	 * 
	 * {@inheritDoc} 
	 */
	@Override
	public void addObservations(double[] source, double[] destination,
			int startTime, int numTimeSteps) throws Exception {
		double[][] sourceMatrix = new double[numTimeSteps][1];
		MatrixUtils.copyIntoColumn(sourceMatrix, 0, 0, source, startTime, numTimeSteps);
		double[][] destMatrix = new double[destination.length][1];
		MatrixUtils.copyIntoColumn(destMatrix, 0, 0, destination, startTime, numTimeSteps);
		addObservations(sourceMatrix, destMatrix);
	}

	@Override
	public void addObservations(double[][] source, double[][] destination) throws Exception {
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		int thisSourceDimensions = source[0].length;
		int thisDestDimensions = destination[0].length;
		if ((thisDestDimensions != destDimensions) || (thisSourceDimensions != sourceDimensions)) {
			throw new Exception("Cannot add observsations for source and destination variables " +
					" of " + thisSourceDimensions + " and " + thisDestDimensions +
					" dimensions respectively for TE calculator set up for " + sourceDimensions + " " +
					destDimensions + " source and destination dimensions respectively");
		}
		if (vectorOfJointSourceObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		vectorOfJointSourceObservations.add(source);
		vectorOfJointDestinationObservations.add(destination);
	}
	
	@Override
	public void addObservations(double[][] source, double[][] destination,
			int startTime, int numTimeSteps) throws Exception {
		double[][] sourceToAdd = new double[numTimeSteps][source[0].length];
		System.arraycopy(source, startTime, sourceToAdd, 0, numTimeSteps);
		double[][] destToAdd = new double[numTimeSteps][destination[0].length];
		System.arraycopy(destination, startTime, destToAdd, 0, numTimeSteps);
		addObservations(sourceToAdd, destToAdd);
	}

	@Override
	public void finaliseAddObservations() {
		// First work out the size to allocate the joint vectors, and do the allocation:
		totalObservations = 0;
		for (double[][] destination : vectorOfJointDestinationObservations) {
			totalObservations += destination.length - k;
		}
		destPastVectors = new double[totalObservations][k * destDimensions];
		destNextVectors = new double[totalObservations][destDimensions];
		sourceVectors = new double[totalObservations][sourceDimensions];
		
		// Construct the joint vectors from the given observations
		int startObservation = 0;
		Iterator<double[][]> iterator = vectorOfJointDestinationObservations.iterator();
		for (double[][] source : vectorOfJointSourceObservations) {
			double[][] destination = iterator.next();
			double[][] currentDestPastVectors = makeJointVectorForPast(destination);
			MatrixUtils.arrayCopy(currentDestPastVectors, 0, 0,
					destPastVectors, startObservation, 0, currentDestPastVectors.length,
					k * destDimensions);
			MatrixUtils.arrayCopy(destination, k, 0,
					destNextVectors, startObservation, 0,
					destination.length - k, destDimensions);
			MatrixUtils.arrayCopy(source, k - 1, 0,
					sourceVectors, startObservation, 0,
					source.length - k, sourceDimensions);
			startObservation += destination.length - k;
		}
		
		// Now set the joint vectors in the kernel estimators
		teKernelEstimator.setObservations(destPastVectors, destNextVectors, sourceVectors);

		// Store whether there was more than one observation set:
		addedMoreThanOneObservationSet = vectorOfJointDestinationObservations.size() > 1;
		
		// And clear the vector of observations
		vectorOfJointSourceObservations = null;
		vectorOfJointDestinationObservations = null;
	}
	
	@Override
	public double computeAverageLocalOfObservations() throws Exception {
		
		double te = 0.0;
		if (debug) {
			MatrixUtils.printMatrix(System.out, destPastVectors);
		}
		for (int b = 0; b < totalObservations; b++) {
			TransferEntropyKernelCounts kernelCounts = teKernelEstimator.getCount(destPastVectors[b],
					destNextVectors[b], sourceVectors[b], b);
			double logTerm = 0.0;
			double cont = 0.0;
			if (kernelCounts.countNextPastSource > 0) {
				logTerm = ((double) kernelCounts.countNextPastSource / (double) kernelCounts.countPastSource) /
							((double) kernelCounts.countNextPast / (double) kernelCounts.countPast);
				cont = Math.log(logTerm);
			}
			te += cont;
			if (debug) {
				System.out.println(b + ": " + destPastVectors[b][0] + " (" +
						kernelCounts.countNextPastSource + " / " + kernelCounts.countPastSource + ") / (" +
						kernelCounts.countNextPast + " / " + kernelCounts.countPast + ") = " + 
						logTerm + " -> " + (cont/log2) + " -> sum: " + (te/log2));
			}
		}
		lastAverage = te / (double) totalObservations / log2;
		return lastAverage;
	}

	/**
	 * <p>Computes the average Transfer Entropy for the previously supplied observations,
	 *  using the Grassberger correction for the point count k: log_e(k) ~= digamma(k).</p>
	 *  
	 * <p><b>HOWEVER</b> -- Kaiser and Schreiber, Physica D 166 (2002) pp. 43-62 suggest (on p. 57) that for the TE
	 *   though the adverse correction of the bias correction is worse than the correction
	 *   itself (because the probabilities being multiplied/divided are not independent), 
	 *   so recommend not to use this method.
	 * </p>
	 * <p>As such, it is implemented here for testing purposes only.</p>
	 * 
	 * @see #computeAverageLocalOfObservations()
	 */
	public double computeAverageLocalOfObservationsWithCorrection() throws Exception {
		
		double te = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			TransferEntropyKernelCounts kernelCounts = teKernelEstimator.getCount(destPastVectors[b],
					destNextVectors[b], sourceVectors[b], b);
			double cont = 0.0;
			if (kernelCounts.countNextPastSource > 0) {
				cont = MathsUtils.digamma(kernelCounts.countNextPastSource) - 
						MathsUtils.digamma(kernelCounts.countPastSource) -
						MathsUtils.digamma(kernelCounts.countNextPast) +
						MathsUtils.digamma(kernelCounts.countPast);
			}
			te += cont;
			/*
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + (cont/log2) + " -> sum: " + (te/log2));
			}
			*/
		}
		// Average it, and convert results to bits
		lastAverage = te / (double) totalObservations / log2;
		return lastAverage;
	}

	@Override
	public double[] computeLocalOfPreviousObservations() throws Exception {
		return computeLocalUsingPreviousObservations(null, null, true);
	}

	/**
	 * {@inheritDoc} 
	 * 
	 * I'm not convinced it is such a good idea to do this 
	 * on data not used for the PDFs with a
	 * kernel estimator (since one can now get kernel estimates
	 * for probabilities of zero now) but I've implemented
	 * it anyway. I guess getting kernel estimates of zero here
	 * is no different than what
	 * can occur with dynamic correlation exclusion.
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] source, double[][] destination) throws Exception {
		return computeLocalUsingPreviousObservations(source, destination, false);
	}
	
	/**
	 * Computes local values for observations of a single-dimensional source and destination pair.
	 * This call is only allowed if the source and destination
	 * dimensions are set to 1.
	 * 
	 * {@inheritDoc} 
	 * 
	 * I'm not convinced it is such a good idea to do this 
	 * on data not used for the PDFs with a
	 * kernel estimator (since one can now get kernel estimates
	 * for probabilities of zero now) but I've implemented
	 * it anyway. I guess getting kernel estimates of zero here
	 * is no different than what
	 * can occur with dynamic correlation exclusion.
	 */
	@Override
	public double[] computeLocalUsingPreviousObservations(
			double[] newSourceObservations, double[] newDestObservations)
			throws Exception {
		double[][] sourceMatrix = new double[newSourceObservations.length][1];
		MatrixUtils.copyIntoColumn(sourceMatrix, 0, newSourceObservations);
		double[][] destMatrix = new double[newDestObservations.length][1];
		MatrixUtils.copyIntoColumn(destMatrix, 0, newDestObservations);
		return computeLocalUsingPreviousObservations(sourceMatrix, destMatrix);
	}
	
	/**
	 * Private utility function to implement
	 * {@link #computeLocalOfPreviousObservations()} and
	 * {@link #computeLocalUsingPreviousObservations(double[][], double[][])}
	 * 
	 * @param source
	 * @param destination
	 * @param isPreviousObservations
	 * @return
	 * @throws Exception
	 */
	private double[] computeLocalUsingPreviousObservations(double[][] source,
			double[][] destination, boolean isPreviousObservations) throws Exception {
		
		double[][] newDestPastVectors;
		double[][] newDestNextValues;
		double[][] newSourceValues;

		if (isPreviousObservations) {
			// We've already computed the joint vectors for these observations
			newDestPastVectors = destPastVectors;
			newDestNextValues = destNextVectors;
			newSourceValues = sourceVectors;
		} else {
			// We need to compute a new set of joint vectors
			newDestPastVectors = makeJointVectorForPast(destination);
			newDestNextValues = new double[destination.length - k][destDimensions];
			MatrixUtils.arrayCopy(destination, k, 0,
					newDestNextValues, 0, 0,
					destination.length - k, destDimensions);
			newSourceValues = new double[source.length - k][sourceDimensions];
			MatrixUtils.arrayCopy(source, k - 1, 0,
					newSourceValues, 0, 0,
					source.length - k, sourceDimensions);
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
		localProbNextCondPast = new double[numLocalObservations];
		double avKernelCount = 0;
		TransferEntropyKernelCounts kernelCounts;
		for (int b = 0; b < numLocalObservations; b++) {
			// System.out.print("Observation number " + String.valueOf(b) + "\n");
			if (isPreviousObservations) {
				kernelCounts = teKernelEstimator.getCount(
						newDestPastVectors[b],
						newDestNextValues[b], newSourceValues[b], b);
			} else {
				kernelCounts = teKernelEstimator.getCount(
						newDestPastVectors[b],
						newDestNextValues[b], newSourceValues[b], -1);
			}
			avKernelCount += kernelCounts.countNextPastSource;
			double logTerm = 0.0;
			double local = 0.0;
			if (kernelCounts.countPast > 0) {
				// Store this ratio for a potential active info calculation later
				localProbNextCondPast[b] = (double) kernelCounts.countNextPast / (double) kernelCounts.countPast;
			}
			if (kernelCounts.countNextPastSource > 0) {
				logTerm = ((double) kernelCounts.countNextPastSource / (double) kernelCounts.countPastSource) /
							localProbNextCondPast[b];
				local = Math.log(logTerm) / log2;
			}
			localTE[offset + b] = local;
			te += local;
			/*
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + (local) + " -> sum: " + (te));
			}
			*/
		}
		avKernelCount = avKernelCount / (double) numLocalObservations;
		if (debug) {
			System.out.printf("Average kernel count was %.3f\n", avKernelCount);
		}
		lastAverage = te / (double) numLocalObservations;
		return localTE;
	}

	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(totalObservations, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}
	
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception {
		
		int numPermutationsToCheck = newOrderings.length;
		
		double actualTE = computeAverageLocalOfObservations();
		
		// Space for the source observations:
		double[][] oldSourceValues = sourceVectors;
		
		// Check that the largest value in the newOrderings is within range:
		int maxNewIndex = MatrixUtils.max(newOrderings);
		if (maxNewIndex >= sourceVectors.length) {
			throw new Exception("Cannot prescribe a new ordering index of " + maxNewIndex +
					" since this is outside the range 0..n-1, where n=" +
					sourceVectors.length + " is the number of observations that " +
					"have been supplied to the calculator.");
		}
		
		int countWhereTeIsMoreSignificantThanOriginal = 0;
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		for (int p = 0; p < numPermutationsToCheck; p++) {
			// Generate a new re-ordered data set for the source in the destPastSourceVectors 
			//  and destNextPastSourceVectors vectors
			sourceVectors = MatrixUtils.extractSelectedTimePoints(oldSourceValues, newOrderings[p]);
			
			// Make the equivalent operations of intialise
			teKernelEstimator.initialise(k * destDimensions,
					sourceDimensions, kernelWidth, kernelWidth);
			// Make the equivalent operations of setObservations:
			teKernelEstimator.setObservations(destPastVectors, destNextVectors, sourceVectors);
			// And get a TE value for this realisation:
			double newTe = computeAverageLocalOfObservations();
			measDistribution.distribution[p] = newTe;
			if (newTe >= actualTE) {
				countWhereTeIsMoreSignificantThanOriginal++;
			}
		}
		
		// Restore the local variables:
		lastAverage = actualTE;
		sourceVectors = oldSourceValues;
		// And set the kernel estimator back to their previous state
		teKernelEstimator.initialise(k * destDimensions,
				sourceDimensions, kernelWidth, kernelWidth);
		teKernelEstimator.setObservations(destPastVectors, destNextVectors, sourceVectors);
		
		// And return the significance
		measDistribution.pValue = (double) countWhereTeIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualTE;
		return measDistribution;
	}

	/**
	 * Computes the local active info storage for the previous supplied observations.
	 * 
	 * <p>Where more than one time series has been added, the array
	 *  contains the local values for each tuple in the order in
	 *  which they were added.</p>
	 * 
	 * <p>If there was only a single time series added, the array
	 *  contains k zero values before the local values.
	 *  (This means the length of the return array is the same
	 *  as the length of the input time series).
	 *  </p>
	 *  
	 * <p>Precondition: The user must have computed the local TEs first for these
	 *  vectors.</p>
	 * 
	 * @see ActiveInfoStorageCalculatorKernel#computeLocalOfPreviousObservations()
	 */
	public double[] computeLocalActiveOfPreviousObservations() throws Exception {
		return computeLocalActiveUsingPreviousObservations(null, true);
	}

	/**
	 * Computes local active info storage for the given observations, using the previously supplied
	 * observations to compute the PDFs.
	 * 
	 * <p>I don't think it's such a good idea to do this for continuous variables (e.g. where
	 * one can get kernel estimates for probabilities of zero now) but I've implemented
	 * it anyway. I guess getting kernel estimates of zero here is no different than what
	 * can occur with dynamic correlation exclusion.</p>
	 * 
	 * <p>Precondition: The user must have computed the local TEs first for these
	 *  vectors</p>
	 * 
	 * @param source
	 * @param destination
	 * @return
	 * @throws Exception
	 * @see ActiveInfoStorageCalculatorKernel#computeLocalUsingPreviousObservations(double[])
	 */
	public double[] computeLocalActiveUsingPreviousObservations(double[][] destination) throws Exception {
		return computeLocalActiveUsingPreviousObservations(destination, false);
	}
	
	/**
	 * Private utility function to implement
	 *  {@link #computeLocalActiveOfPreviousObservations()} and
	 *  {@link #computeLocalActiveUsingPreviousObservations(double[][])}
	 * 
	 * @param source
	 * @param destination
	 * @param isPreviousObservations
	 * @return
	 * @throws Exception
	 */
	private double[] computeLocalActiveUsingPreviousObservations(
			double[][] destination, boolean isPreviousObservations) throws Exception {

		// Precondition: the local TE must have already been computed
		if (localProbNextCondPast == null) {
			throw new RuntimeException("A local TE must have been computed before " +
					"the local active info storage can be computed by TransferEntropyCalculatorMultiVariateKernel");
		}
		
		double[][] newDestNextValues;

		// Set the observations on the kernel estimator:
		nextStateKernelEstimator.setObservations(destNextVectors);
		
		// Now set which observations we're going to compute the local 
		//  active info of:
		if (isPreviousObservations) {
			// We've already computed the joint vectors for these observations
			newDestNextValues = destNextVectors;
		} else {
			// We need to compute a new set of joint vectors
			newDestNextValues = new double[destination.length - k][destDimensions];
			MatrixUtils.arrayCopy(destination, k, 0,
					newDestNextValues, 0, 0,
					destination.length - k, destDimensions);
		}

		int numLocalObservations = newDestNextValues.length;
		double[] localActive;
		int offset = 0;
		if (isPreviousObservations && addedMoreThanOneObservationSet) {
			// We're returning the local values for a set of disjoint
			//  observations. So we don't add k zeros to the start
			localActive = new double[numLocalObservations];
			offset = 0;
		} else {
			localActive = new double[numLocalObservations + k];
			offset = k;
		}
		double nextStateProb;
		for (int b = 0; b < numLocalObservations; b++) {
			if (isPreviousObservations) {
				nextStateProb = nextStateKernelEstimator.getProbability(
						newDestNextValues[b], b);
			} else {
				nextStateProb = nextStateKernelEstimator.getProbability(
						newDestNextValues[b], -1);
			}
			double logTerm = 0.0;
			double local = 0.0;
			if (localProbNextCondPast[b] > 0) {
				logTerm = localProbNextCondPast[b] / nextStateProb;
				local = Math.log(logTerm);
			}
			localActive[offset + b] = local / log2;
			/*
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + (local/log2) + " -> sum: " + (te/log2));
			}
			*/
		}
		return localActive;
	}
	
	@Override
	public void setDebug(boolean debug) {
		super.setDebug(debug);
		teKernelEstimator.setDebug(debug);
	}

	/**
	 * Generate an embedding vector for each time step, containing the past k states of the destination.
	 * Note that each state of the destination is a joint vector of destDimensions variables.
	 * Does not include a vector for the first k time steps.
	 * 
	 * @param destination destination time-series
	 * @return array of embedding vectors for each time step of
	 * 	length destination.length - k
	 */
	private double[][] makeJointVectorForPast(double[][] destination) {
		try {
			// We want one less delay vector here - we don't need the last k point,
			//  because there is no next state for these.
			return MatrixUtils.makeDelayEmbeddingVector(destination, k, k-1, destination.length - k);
		} catch (Exception e) {
			// The parameters for the above call should be fine, so we don't expect to
			//  throw an Exception here - embed in a RuntimeException if it occurs 
			throw new RuntimeException(e);
		}
	}
}
