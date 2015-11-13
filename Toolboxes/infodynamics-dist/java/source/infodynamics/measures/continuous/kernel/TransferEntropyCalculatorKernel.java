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
import infodynamics.measures.continuous.kernel.TransferEntropyKernelCounts;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

import java.util.Iterator;

/**
 * <p>Computes the differential transfer entropy (TE) between two univariate
 *  <code>double[]</code> time-series of observations
 *  using box-kernel estimation.
 *  For details on box-kernel estimation, see Kantz and Schreiber (below).</p>
 *  
 * <p>TE was defined by Schreiber (below).
 *  This estimator is realised here by extending
 *  {@link TransferEntropyCommon}.</p>
 *  
 * <p>This is our <b>main</b> class for TE by box-kernel estimation, implementing
 * dynamic correlation exclusion, and optimisation for the box counting.
 * Note the existence of several other classes using box-kernel estimation,
 * which implement TE in various less efficient ways, e.g. 
 * {@link TransferEntropyCalculatorKernelPlain}, {@link TransferEntropyCalculatorKernelPlainIterators}</p>
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
 *  <li>H. Kantz and T. Schreiber, "Nonlinear Time Series Analysis"
 *  (Cambridge University Press, Cambridge, MA, 1997).</li>
 * </ul>
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class TransferEntropyCalculatorKernel
	extends TransferEntropyCommon implements TransferEntropyCalculator {

	protected KernelEstimatorTransferEntropy teKernelEstimator = null;

	// Keep joint vectors so we don't need to regenerate them
	protected double[][] destPastVectors;
	protected double[]   destNextValues;
	protected double[]   sourceValues;
	
	private boolean normalise = true;
	/**
	 * Property name for whether to normalise the incoming variables 
	 * to mean 0, standard deviation 1, or not (default false)
	 */
	public static final String NORMALISE_PROP_NAME = "NORMALISE";
	
	private boolean dynCorrExcl = false;
	private int dynCorrExclTime = 0;
	/**
	 * Property name for a dynamics exclusion time window (see Kantz and Schreiber),
	 * default is 0 which means no dynamic exclusion window.
	 */
	public static final String DYN_CORR_EXCL_TIME_NAME = "DYN_CORR_EXCL";
	
	private boolean forceCompareToAll = false;
	/**
	 * Property name for whether to force the underlying kernel estimators to compare
	 *  each data point to each other (or else allow it to use optimisations)
	 */
	public static final String FORCE_KERNEL_COMPARE_TO_ALL = "FORCE_KERNEL_COMPARE_TO_ALL";
	
	/**
	 * Default value for kernel width
	 */
	public static final double DEFAULT_KERNEL_WIDTH = 0.25;
	/**
	 * Kernel width
	 */
	private double kernelWidth = DEFAULT_KERNEL_WIDTH;
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
	public TransferEntropyCalculatorKernel() {
		super();
		teKernelEstimator = new KernelEstimatorTransferEntropy();
		teKernelEstimator.setNormalise(normalise);
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
		super.initialise(k); // calls initialise();
	}

	@Override
	public void initialise() {
		teKernelEstimator.initialise(k, kernelWidth);
		destPastVectors = null;
		destNextValues = null;
		sourceValues = null;
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
	 * 		<li>{@link #DYN_CORR_EXCL_TIME_NAME} -- a dynamics exclusion time window,
	 * 			also known as Theiler window (see Kantz and Schreiber);
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
		} else if (propertyName.equalsIgnoreCase(DYN_CORR_EXCL_TIME_NAME)) {
			dynCorrExclTime = Integer.parseInt(propertyValue);
			dynCorrExcl = (dynCorrExclTime > 0);
			if (dynCorrExcl) {
				teKernelEstimator.setDynamicCorrelationExclusion(dynCorrExclTime);
			} else {
				teKernelEstimator.clearDynamicCorrelationExclusion();
			}
		} else if (propertyName.equalsIgnoreCase(FORCE_KERNEL_COMPARE_TO_ALL)) {
			forceCompareToAll = Boolean.parseBoolean(propertyValue);
			teKernelEstimator.setForceCompareToAll(forceCompareToAll);
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
	public String getProperty(String propertyName) {
		if (propertyName.equalsIgnoreCase(KERNEL_WIDTH_PROP_NAME) ||
				propertyName.equalsIgnoreCase(EPSILON_PROP_NAME)) {
			return Double.toString(kernelWidth);
		} else if (propertyName.equalsIgnoreCase(NORMALISE_PROP_NAME)) {
			return Boolean.toString(normalise);
		} else if (propertyName.equalsIgnoreCase(DYN_CORR_EXCL_TIME_NAME)) {
			return Integer.toString(dynCorrExclTime);
		} else if (propertyName.equalsIgnoreCase(FORCE_KERNEL_COMPARE_TO_ALL)) {
			return Boolean.toString(forceCompareToAll);
		} else {
			// try the superclass:
			return super.getProperty(propertyName);
		}
	}

	@Override
	public void finaliseAddObservations() {
		// First work out the size to allocate the joint vectors, and do the allocation:
		totalObservations = 0;
		for (double[] destination : vectorOfDestinationObservations) {
			totalObservations += destination.length - k;
		}
		destPastVectors = new double[totalObservations][k];
		destNextValues = new double[totalObservations];
		sourceValues = new double[totalObservations];
		
		// Construct the joint vectors from the given observations
		int startObservation = 0;
		Iterator<double[]> iterator = vectorOfDestinationObservations.iterator();
		for (double[] source : vectorOfSourceObservations) {
			double[] destination = iterator.next();
			double[][] currentDestPastVectors = makeJointVectorForPast(destination);
			MatrixUtils.arrayCopy(currentDestPastVectors, 0, 0,
					destPastVectors, startObservation, 0, currentDestPastVectors.length, k);
			System.arraycopy(destination, k, destNextValues, startObservation, destination.length - k);
			System.arraycopy(source, k - 1, sourceValues, startObservation, source.length - k);
			startObservation += destination.length - k;
		}
		
		// Now set the joint vectors in the kernel estimators
		teKernelEstimator.setObservations(destPastVectors, destNextValues, sourceValues);

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
	
	@Override
	public double computeAverageLocalOfObservations() throws Exception {
		
		double te = 0.0;
		if (debug) {
			MatrixUtils.printMatrix(System.out, destPastVectors);
		}
		for (int b = 0; b < totalObservations; b++) {
			TransferEntropyKernelCounts kernelCounts = teKernelEstimator.getCount(destPastVectors[b],
					destNextValues[b], sourceValues[b], b);
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
						logTerm + " -> " + (cont/Math.log(2.0)) + " -> sum: " + (te/Math.log(2.0)));
			}
		}
		lastAverage = te / (double) totalObservations / Math.log(2.0);
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
	 *   (Bias correction is implemented properly in the Kraskov et al. 
	 *   estimators, see {@link infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov}
	 * </p>
	 * <p>As such, it is implemented here for testing purposes only.</p>
	 * 
	 * @see #computeAverageLocalOfObservations()
	 */
	public double computeAverageLocalOfObservationsWithCorrection() throws Exception {
		
		double te = 0.0;
		int numNoNeighbours = 0;
		double contributionsNoNeighbours = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			TransferEntropyKernelCounts kernelCounts = teKernelEstimator.getCount(destPastVectors[b],
					destNextValues[b], sourceValues[b], b);
			double cont = 0.0;
			// Original code:
			/* if (kernelCounts.countNextPastSource > 0) {
				cont = MathsUtils.digamma(kernelCounts.countNextPastSource) - 
						MathsUtils.digamma(kernelCounts.countPastSource) -
						MathsUtils.digamma(kernelCounts.countNextPast) +
						MathsUtils.digamma(kernelCounts.countPast);
			} */
			// But Schreiber confirmed to me that with dynamic correlation 
			//  exclusion, when you may have no nearest neighbours,
			//  he was allowing contributions from other groups (e.g. countPastSource)
			//  to be added even if the full joint count was zero.
			// Implement it like this:
			//  TODO Do we need to correct the totalObservations
			//   divisor for each digamma sum to account for this?
			//   (Schreiber does not do that; for the moment we'll accept that
			//    as the right approach)
			if (kernelCounts.countPastSource > 0) {
				cont -= MathsUtils.digamma(kernelCounts.countPastSource);
			}
			if (kernelCounts.countNextPast > 0) {
				cont -= MathsUtils.digamma(kernelCounts.countNextPast);
			}
			if (kernelCounts.countPast > 0) {
				cont += MathsUtils.digamma(kernelCounts.countPast);
			}
			if (kernelCounts.countNextPastSource > 0) {
				cont += MathsUtils.digamma(kernelCounts.countNextPastSource);
			} else {
				// These contributions are from a set with no neighbours in 
				//  the full joint space.
				contributionsNoNeighbours += cont;
				numNoNeighbours++;
			}
			te += cont;
			/*
			if (debug) {
				System.out.println(b + ": " + cont + " -> " + (cont/Math.log(2.0)) + " -> sum: " + (te/Math.log(2.0)));
			}
			*/
		}
		// Average it, and convert results to bytes
		lastAverage = te / (double) totalObservations / Math.log(2.0);
		if (debug) {
			System.out.printf("TE=%.4f, with %d contributions from 0 neighbour sets being %.4f\n",
					lastAverage,
					numNoNeighbours,
					contributionsNoNeighbours / (double) totalObservations / Math.log(2.0));
		}
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
	public double[] computeLocalUsingPreviousObservations(double[] source, double[] destination) throws Exception {
		return computeLocalUsingPreviousObservations(source, destination, false);
	}
	
	/**
	 * Private utility function to implement
	 * {@link #computeLocalOfPreviousObservations()} and
	 * {@link #computeLocalUsingPreviousObservations(double[], double[])}
	 * 
	 * @param source
	 * @param destination
	 * @param isPreviousObservations
	 * @return
	 * @throws Exception
	 */
	private double[] computeLocalUsingPreviousObservations(double[] source, double[] destination, boolean isPreviousObservations) throws Exception {
		
		double[][] newDestPastVectors;
		double[] newDestNextValues;
		double[] newSourceValues;

		if (isPreviousObservations) {
			// We've already computed the joint vectors for these observations
			newDestPastVectors = destPastVectors;
			newDestNextValues = destNextValues;
			newSourceValues = sourceValues;
		} else {
			// We need to compute a new set of joint vectors
			newDestPastVectors = makeJointVectorForPast(destination);
			newDestNextValues = MatrixUtils.select(destination, k, destination.length - k);
			newSourceValues = MatrixUtils.select(source, k - 1, source.length - k);
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
		TransferEntropyKernelCounts kernelCounts;
		for (int b = 0; b < numLocalObservations; b++) {
			if (isPreviousObservations) {
				kernelCounts = teKernelEstimator.getCount(newDestPastVectors[b],
						newDestNextValues[b], newSourceValues[b], b);
			} else {
				kernelCounts = teKernelEstimator.getCount(newDestPastVectors[b],
						newDestNextValues[b], newSourceValues[b], -1);
			}
			double logTerm = 0.0;
			double local = 0.0;
			if (kernelCounts.countNextPastSource > 0) {
				logTerm = ((double) kernelCounts.countNextPastSource / (double) kernelCounts.countPastSource) /
							((double) kernelCounts.countNextPast / (double) kernelCounts.countPast);
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
		double[] oldSourceValues = sourceValues;
		
		int countWhereTeIsMoreSignificantThanOriginal = 0;
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		for (int p = 0; p < numPermutationsToCheck; p++) {
			// Generate a new re-ordered data set for the source in the destPastSourceVectors 
			//  and destNextPastSourceVectors vectors
			sourceValues = MatrixUtils.extractSelectedTimePoints(oldSourceValues, newOrderings[p]);
			
			// Make the equivalent operations of intialise
			teKernelEstimator.initialise(k, kernelWidth);
			// Make the equivalent operations of setObservations:
			teKernelEstimator.setObservations(destPastVectors, destNextValues, sourceValues);
			// And get a TE value for this realisation:
			double newTe = computeAverageLocalOfObservations();
			measDistribution.distribution[p] = newTe;
			if (newTe >= actualTE) {
				countWhereTeIsMoreSignificantThanOriginal++;
			}
		}
		
		// Restore the local variables:
		lastAverage = actualTE;
		sourceValues = oldSourceValues;
		// And set the kernel estimator back to their previous state
		teKernelEstimator.initialise(k, kernelWidth);
		teKernelEstimator.setObservations(destPastVectors, destNextValues, sourceValues);
		
		// And return the significance
		measDistribution.pValue = (double) countWhereTeIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualTE;
		return measDistribution;
	}

	@Override
	public void setDebug(boolean debug) {
		super.setDebug(debug);
		teKernelEstimator.setDebug(debug);
	}
}
