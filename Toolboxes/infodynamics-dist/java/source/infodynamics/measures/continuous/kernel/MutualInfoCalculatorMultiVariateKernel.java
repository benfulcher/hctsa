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

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.MutualInfoMultiVariateCommon;

/**
 * <p>Computes the differential mutual information of two given multivariate sets of
 *  observations (implementing {@link MutualInfoCalculatorMultiVariate}),
 *  using box-kernel estimation.
 *  For details on box-kernel estimation, see Kantz and Schreiber (below).</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link MutualInfoCalculatorMultiVariate},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #MutualInfoCalculatorMultiVariateKernel()}.</li>
 *  <li>Further properties are available, see {@link #setProperty(String, String)};</li>
 *  <li>An additional {@link #initialise(int, int, double)} option;</li>
 *  <li>Additional utility methods for computing other information-theoretic values
 *  are available here (e.g. {@link #computeAverageJointEntropy()}) which
 *  can be called after all observations are supplied.</li>
 *  </ul>
 * </p>
 * 
 * <p>
 * TODO Use only a single kernel estimator class for the joint space, and compute other
 *  probabilities from this. This will save much time.
 * </p>
 * 
 * @see "H. Kantz and T. Schreiber, 'Nonlinear Time Series Analysis'.
 *   Cambridge, MA: Cambridge University Press, 1997"
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class MutualInfoCalculatorMultiVariateKernel
	extends MutualInfoMultiVariateCommon
	implements MutualInfoCalculatorMultiVariate, Cloneable {

	protected KernelEstimatorMultiVariate mvkeSource = null;
	protected KernelEstimatorMultiVariate mvkeDest = null;
	protected KernelEstimatorMultiVariate mvkeJoint = null;

	private boolean normalise = true;
	/**
	 * Property name for whether to normalise the incoming variables 
	 * to mean 0, standard deviation 1, or not (default false)
	 */
	public static final String NORMALISE_PROP_NAME = "NORMALISE";
	
	private boolean dynCorrExcl = false;
	private int dynCorrExclTime = 100;
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
	 * Kernel width currently in use
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
	 * Construct an instance of the box-kernel MI calculator
	 */
	public MutualInfoCalculatorMultiVariateKernel() {
		// Create our kernel estimator objects:
		mvkeSource = new KernelEstimatorMultiVariate();
		mvkeDest = new KernelEstimatorMultiVariate();
		mvkeJoint = new KernelEstimatorMultiVariate();
		mvkeSource.setNormalise(normalise);
		mvkeDest.setNormalise(normalise);
		mvkeJoint.setNormalise(normalise);
	}

	public void initialise(int sourceDimensions, int destDimensions) {
		initialise(sourceDimensions, destDimensions, kernelWidth);
	}

	/**
	 * Initialise the calculator for (re-)use, with a specific kernel width
	 * with number of source
	 * and destination joint variables specified, and existing
	 * (or default) values of other parameters,.
	 * Clears an PDFs of previously supplied observations.
	 *
	 * @param sourceDimensions the number of joint variables in the source
	 * @param destDimensions the number of joint variables in the destination
	 * @param kernelWidth if {@link #NORMALISE_PROP_NAME} property has
	 *  been set, then this kernel width corresponds to the number of
	 *  standard deviations from the mean (otherwise it is an absolute value)
	 */
	public void initialise(int sourceDimensions, int destDimensions, double kernelWidth) {
		super.initialise(sourceDimensions, destDimensions);
		// Store kernel width for local use
		this.kernelWidth = kernelWidth;
		// Initialise the kernel estimators:
		mvkeSource.initialise(sourceDimensions, kernelWidth);
		mvkeDest.initialise(destDimensions, kernelWidth);
		mvkeJoint.initialise(sourceDimensions + destDimensions, kernelWidth);
	}

	public void finaliseAddObservations() {
		// Get the observations properly stored in the sourceObservations[][] and
		//  destObservations[][] arrays.
		try {
			// Currently, the throws declaration in super is only there to
			//  allow other children to throw Exceptions - there are no compile
			//  time Exceptions being thrown there.
			super.finaliseAddObservations();
		} catch (Exception e) {
			// So we cast any found Exception to a RuntimeException
			throw new RuntimeException(e);
		}

		// Now assign these observations to the underlying kernel estimators:
		mvkeSource.setObservations(sourceObservations);
		mvkeDest.setObservations(destObservations);
		try {
			// This will only throw an exception, in theory, if
			//  the time length of the observations is different
			//  or if they are null - since we constructed them 
			//  neither of these should be the case.
			mvkeJoint.setObservations(sourceObservations, destObservations);
		} catch (Exception e) {
			throw new RuntimeException("Unhandled exception from MultivariateKernelEstimator.setObservations(double[][], double[][])", e);
		}
	}

	/**
	 * Compute the average MI from the previously supplied observations.
	 * 
	 * @return the average MI in bits
	 */
	public double computeAverageLocalOfObservations() {
		double mi = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob1 = mvkeSource.getProbability(sourceObservations[b], b);
			double prob2 = mvkeDest.getProbability(destObservations[b], b);
			double probJoint = mvkeJoint.getProbability(sourceObservations[b], destObservations[b], b);
			double logTerm = 0.0;
			double cont = 0.0;
			if (probJoint > 0.0) {
				// If we have counted joint correlations, we must have marginals for each
				logTerm = probJoint / (prob1 * prob2);
				cont = Math.log(logTerm);
			}
			mi += cont;
			if (debug) {
				System.out.printf("%d: (%.5f, %.5f, %.5f) %.5f -> %.5f -> %.5f\n",
						b, prob1, prob2, probJoint, logTerm, cont, mi);
			}
		}
		lastAverage = mi / (double) totalObservations / Math.log(2.0);
		miComputed = true;
		return lastAverage;
	}

	/**
	 * Extra utility method to return the joint entropy, for the source
	 * and destination variables considered jointly, using the previously supplied
	 * observations.
	 * 
	 * @return the average joint entropy in bits.
	 */
	public double computeAverageJointEntropy() {
		double entropy = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob = mvkeJoint.getProbability(sourceObservations[b], destObservations[b], b);
			double cont = 0.0;
			if (prob > 0.0) {
				cont = - Math.log(prob);
			}
			entropy += cont;
			if (debug) {
				System.out.println(b + ": " + prob + " -> " + cont/Math.log(2.0) + " -> sum: " + (entropy/Math.log(2.0)));
			}
		}
		return entropy / (double) totalObservations / Math.log(2.0);
	}

	/**
	 * Extra utility method to return the entropy of the (first) source variable
	 * 
	 * @return entropy of source variable in bits.
	 */
	public double computeAverageEntropyOfObservation1() {
		double entropy = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob = mvkeSource.getProbability(sourceObservations[b], b);
			double cont = 0.0;
			// Comparing the prob to 0.0 should be fine - it would have to be
			//  an impossible number of samples for us to hit machine resolution here.
			if (prob > 0.0) {
				cont = -Math.log(prob);
			}
			entropy += cont;
			if (debug) {
				System.out.println(b + ": " + prob + " -> " + cont/Math.log(2.0) + " -> sum: " + (entropy/Math.log(2.0)));
			}
		}
		return entropy / (double) totalObservations / Math.log(2.0);
	}

	/**
	 * Extra utility method to return the entropy of the second (destination) variable
	 * 
	 * @return entropy of the destination in bits
	 */
	public double computeAverageEntropyOfObservation2() {
		double entropy = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob = mvkeDest.getProbability(destObservations[b], b);
			double cont = 0.0;
			if (prob > 0.0) {
				cont = -Math.log(prob);
			}
			entropy += cont;
			if (debug) {
				System.out.println(b + ": " + prob + " -> " + cont/Math.log(2.0) + " -> sum: " + (entropy/Math.log(2.0)));
			}
		}
		return entropy / (double) totalObservations / Math.log(2.0);
	}

	/**
	 * Extra utility method to return the information distance between the source 
	 * and destination
	 * 
	 * @return information distance between the source and destination in bits.
	 */
	public double computeAverageInfoDistanceOfObservations() {
		double infoDistance = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob1 = mvkeSource.getProbability(sourceObservations[b], b);
			double prob2 = mvkeDest.getProbability(destObservations[b], b);
			double probJoint = mvkeJoint.getProbability(sourceObservations[b], destObservations[b], b);
			double logTerm = 0.0;
			double cont = 0.0;
			if (probJoint > 0.0) {
				logTerm = (prob1 * prob2) / (probJoint * probJoint);
				cont = Math.log(logTerm);
			}
			infoDistance += cont;
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + (cont/Math.log(2.0)) + " -> sum: " + (infoDistance/Math.log(2.0)));
			}
		}
		return infoDistance / (double) totalObservations / Math.log(2.0);
	}

	/**
	 * <p>Computes the local values of the MI,
	 *  for each valid observation in the previously supplied observations
	 *  (with PDFs computed using all of the previously supplied observation sets).</p>
	 *  
	 * <p>If the samples were supplied via a single call such as
	 * {@link #setObservations(double[])},
	 * then the return value is a single time-series of local
	 * channel measure values corresponding to these samples.</p>
	 * 
	 * <p>Otherwise where disjoint time-series observations were supplied using several 
	 *  calls such as {@link addObservations(double[])}
	 *  then the local values for each disjoint observation set will be appended here
	 *  to create a single "time-series" return array.</p>
	 * 
	 * @return the "time-series" of local MIs in bits
	 * @throws Exception
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		return computeLocalUsingPreviousObservations(sourceObservations, destObservations, true);
	}


	/**
	 * Compute the local MI values for each of the
	 * supplied samples in <code>states1</code> and <code>states2</code>.
	 * 
	 * <p>PDFs are computed using all of the previously supplied
	 * observations, but not those in <code>states1</code> and <code>states2</code>
	 * (unless they were
	 * some of the previously supplied samples).</p>
	 * 
	 * <p>Calls to this method will not harness dynamic correlation exclusion (if set)
	 *  since we don't know whether it's the same time set or not.
	 * </p>
	 *  
	 * @param states1 series of multivariate observations for the source variable
	 *  (first index is time or observation index, second is variable number)
	 * @param states2 series of multivariate observations for the destination variable
	 *  (first index is time or observation index, second is variable number).
	 *  Length must match <code>source</code>, and their indices must correspond.
	 * @return the local values in bits.
	 *  If the {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF}
	 *  property was set to say k, then the local values align with the
	 *  destination value (i.e. after the given delay k). As such, the
	 *  first k values of the array will be zeros. 
	 */
	public double[] computeLocalUsingPreviousObservations(double states1[][], double states2[][]) {
		return computeLocalUsingPreviousObservations(states1, states2, false);
	}
	
	/**
	 * Protected utility function to compute the local MI values for each of the
	 * supplied samples in <code>states1</code> and <code>states2</code>.
	 * 
	 * <p>PDFs are computed using all of the previously supplied
	 * observations. <code>isOurPreviousObservations</code> indicates whether
	 * those in <code>states1</code> and <code>states2</code>
	 * were some of the previously supplied samples.</p>
	 * 
	 * @param states1 provided source observations
	 * @param states2 provided destination observations
	 * @param isOurPreviousObservations whether these are our previous
	 *  observations - this determines whether to add zeros for the first
	 *  timeDiff local values, and also 
	 *  whether to set the internal lastAverage field,
	 *  which is returned by later calls to {@link #getLastAverage()}
	 * @return the local values in bits.
	 *  If the {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF}
	 *  property was set to say k, then the local values align with the
	 *  destination value (i.e. after the given delay k). As such, the
	 *  first k values of the array will be zeros.
	 */
	protected double[] computeLocalUsingPreviousObservations(double states1[][],
			double states2[][], boolean isOurPreviousObservations) {
		
		double mi = 0.0;
		int timeSteps = states1.length;
		double[] localMi = new double[timeSteps];
		double prob1, prob2, probJoint;
		for (int b = 0; b < timeSteps; b++) {
			if (isOurPreviousObservations) {
				// We've been called with our previous observations, so we
				//  can pass the time step through for dynamic correlation exclusion
				prob1 = mvkeSource.getProbability(states1[b], b);
				prob2 = mvkeDest.getProbability(states2[b], b);
				probJoint = mvkeJoint.getProbability(states1[b], states2[b], b);
			} else {
				// We don't know whether these were our previous observation or not
				//  so we don't do dynamic correlation exclusion
				prob1 = mvkeSource.getProbability(states1[b]);
				prob2 = mvkeDest.getProbability(states2[b]);
				probJoint = mvkeJoint.getProbability(states1[b], states2[b]);
			}
			double logTerm = 0.0;
			localMi[b] = 0.0;
			if (probJoint > 0.0) {
				// By necessity prob1 and prob2 will be > 0.0
				logTerm = probJoint / (prob1 * prob2);
				localMi[b] = Math.log(logTerm) / Math.log(2.0);
			}
			mi += localMi[b];
			if (debug) {
				System.out.printf("%d: (%.5f, %.5f, %.5f) %.5f -> %.5f -> %.5f\n",
						b, prob1, prob2, probJoint, logTerm, localMi[b], mi);
			}
		}
		lastAverage = mi / (double) totalObservations;
		miComputed = true;
		return localMi;
	}

	
	/**
	 * Compute the local joint entropy values of the previously provided
	 *  observations.
	 *  
	 * @param states1 provided source observations
	 * @param states2 provided destination observations
	 * @return the local joint entropies in bits
	 */
	public double[] computeLocalJointEntropyOfPreviousObservations() throws Exception {
		return computeLocalJointEntropyUsingPreviousObservations(sourceObservations,
				destObservations, true);
	}

	/**
	 * Compute the local joint entropy values for these given values, using the previously provided
	 *  observations to compute the probabilities.
	 * Calls to this method will not harness dynamic correlation exclusion (if set)
	 *  since we don't know whether it's the same time set or not.
	 *  
	 * @param states1 provided source observations
	 * @param states2 provided destination observations
	 * @return the local joint entropies in bits
	 */
	public double[] computeLocalJointEntropyUsingPreviousObservations(double states1[][], double states2[][]) {
		return computeLocalJointEntropyUsingPreviousObservations(states1,
				states2, false);
	}
	
	/**
	 * Internal implementation
	 * 
	 * @param states1
	 * @param states2
	 * @param isOurPreviousObservations
	 * @return
	 */
	private double[] computeLocalJointEntropyUsingPreviousObservations(
			double states1[][], double states2[][], boolean isOurPreviousObservations) {
		int timeSteps = states1.length;
		double[] localJoint = new double[timeSteps];
		double prob;
		for (int b = 0; b < totalObservations; b++) {
			if (isOurPreviousObservations) {
				prob = mvkeJoint.getProbability(sourceObservations[b], destObservations[b], b);
			} else {
				prob = mvkeJoint.getProbability(sourceObservations[b], destObservations[b]);
			}
			localJoint[b] = 0.0;
			if (prob > 0.0) {
				localJoint[b] = - Math.log(prob) / Math.log(2.0);
			}
			if (debug) {
				System.out.println(b + ": " + prob + " -> " + localJoint[b]);
			}
		}
		return localJoint;
	}

	/**
	 * Compute the local entropy values for the previously provided
	 *  observations for the (first) source variable 
	 *  (using those previous observations to compute the PDFs).
	 * 
	 * @return array of local entropies for the source observations
	 */
	public double[] computeLocalEntropy1OfPreviousObservations() {
		return computeLocalEntropyFromPreviousObservations(sourceObservations, 1, true);
	}

	/**
	 * Compute the local entropy values for these given source values,
	 * using the previously provided
	 *  observations for the (first) source variable to compute the probabilities.
	 * Calls to this method will not harness dynamic correlation exclusion (if set)
	 *  since we don't know whether it's the same time set or not.
	 * 
	 * @param states1 provided source observations
	 *  
	 * @return array of local entropies for these source observations
	 */
	public double[] computeLocalEntropy1UsingPreviousObservations(double[][] states) {
		return computeLocalEntropyFromPreviousObservations(states, 1, false);
	}

	/**
	 * Compute the local entropy values for the previously provided
	 *  observations for the (second) destination variable 
	 *  (using those previous observations to compute the PDFs).
	 * 
	 * @return array of local entropies for the destination observations
	 */
	public double[] computeLocalEntropy2OfPreviousObservations() {
		return computeLocalEntropyFromPreviousObservations(destObservations, 2, true);
	}

	/**
	 * Compute the local entropy values for these given destination values,
	 * using the previously provided
	 *  observations for the (second) destination variable to compute the probabilities.
	 * Calls to this method will not harness dynamic correlation exclusion (if set)
	 *  since we don't know whether it's the same time set or not.
	 * 
	 * @param states2 provided destination observations
	 *  
	 * @return array of local entropies for these destination observations
	 */
	public double[] computeLocalEntropy2UsingPreviousObservations(double states[][]) {
		return computeLocalEntropyFromPreviousObservations(states, 2, false);
	}

	/**
	 * Private utility function to implement {@link #computeLocalEntropy1OfPreviousObservations()},
	 * {@link #computeLocalEntropy1UsingPreviousObservations(double[][])},
	 * {@link #computeLocalEntropy2OfPreviousObservations()} and
	 * {@link #computeLocalEntropy2UsingPreviousObservations(double[][])}.
	 *  
	 * @param states provided set of observations
	 * @param useProbsForWhichVar use 1 for variable 1, 2 for variable 2
	 * @param isOurPreviousObservations whether these are the previously provided
	 * observations or not.
	 * @return array of local entropies for these observations
	 */
	private double[] computeLocalEntropyFromPreviousObservations(
			double states[][], int useProbsForWhichVar, boolean isOurPreviousObservations) {
		int timeSteps = states.length;
		double[] localEntropy = new double[timeSteps];
		double prob;
		for (int b = 0; b < totalObservations; b++) {
			if (useProbsForWhichVar == 1) {
				if (isOurPreviousObservations) {
					prob = mvkeSource.getProbability(states[b], b);
				} else {
					prob = mvkeSource.getProbability(states[b]);
				}
			} else {
				if (isOurPreviousObservations) {
					prob = mvkeDest.getProbability(states[b], b);
				} else {
					prob = mvkeDest.getProbability(states[b]);
				}
			}
			localEntropy[b] = 0.0;
			if (prob > 0.0) {
				localEntropy[b] = - Math.log(prob) / Math.log(2.0);
			}
			if (debug) {
				System.out.println(b + ": " + prob + " -> " + localEntropy[b]);
			}
		}
		return localEntropy;
	}

	/**
	 * Compute the local Info distance values for the previously provided
	 *  observations to compute the probabilities.
	 * 
	 * @return array of local information distances
	 */
	public double[] computeLocalInfoDistanceOfPreviousObservations() {
		return computeLocalInfoDistanceUsingPreviousObservations(sourceObservations,
				destObservations, true);
	}

	/**
	 * Compute the local Info distance values for these given values, using the previously provided
	 *  observations to compute the probabilities.
	 * Calls to this method will not harness dynamic correlation exclusion (if set)
	 *  since we don't know whether it's the same time set or not.
	 * 
	 * @return array of local information distances.
	 */
	public double[] computeLocalInfoDistanceUsingPreviousObservations(double[][] states1, double[][] states2) {
		return computeLocalInfoDistanceUsingPreviousObservations(states1,
				states2, false);
	}

	/**
	 * Protected utility function to implement {@link #computeLocalInfoDistanceOfPreviousObservations()}
	 * and {@link #computeLocalInfoDistanceUsingPreviousObservations(double[][], double[][])}.
	 * 
	 * @param states1 provided source observations
	 * @param states2 provided destination observations
	 * @param isOurPreviousObservations whether these are the previously provided
	 * observations or not.
	 * @return array of local information distances.
	 */
	protected double[] computeLocalInfoDistanceUsingPreviousObservations(
			double[][] states1, double[][] states2, boolean isOurPreviousObservations) {
		int timeSteps = states1.length;
		double[] localInfoDistance = new double[timeSteps];
		double prob1, prob2, probJoint;
		for (int b = 0; b < timeSteps; b++) {
			if (isOurPreviousObservations) {
				prob1 = mvkeSource.getProbability(states1[b], b);
				prob2 = mvkeDest.getProbability(states2[b], b);
				probJoint = mvkeJoint.getProbability(states1[b], states2[b], b);
			} else {
				prob1 = mvkeSource.getProbability(states1[b]);
				prob2 = mvkeDest.getProbability(states2[b]);
				probJoint = mvkeJoint.getProbability(states1[b], states2[b]);
			}
			double logTerm = 0.0;
			localInfoDistance[b] = 0.0;
			if (probJoint > 0.0) {
				logTerm = (prob1 * prob2) / (probJoint * probJoint);
				localInfoDistance[b] = Math.log(logTerm) / Math.log(2.0);
			}
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + localInfoDistance[b]);
			}
		}
		return localInfoDistance;
	}

	/**
	 * <p>Set properties for the kernel MI calculator.
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
	 *  	<li>any valid properties for {@link MutualInfoMultiVariateCommon#setProperty(String, String)}.</li>
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
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(KERNEL_WIDTH_PROP_NAME) ||
				propertyName.equalsIgnoreCase(EPSILON_PROP_NAME)) {
			kernelWidth = Double.parseDouble(propertyValue);
		} else if (propertyName.equalsIgnoreCase(NORMALISE_PROP_NAME)) {
			normalise = Boolean.parseBoolean(propertyValue);
			mvkeSource.setNormalise(normalise);
			mvkeDest.setNormalise(normalise);
			mvkeJoint.setNormalise(normalise);
		} else if (propertyName.equalsIgnoreCase(DYN_CORR_EXCL_TIME_NAME)) {
			dynCorrExclTime = Integer.parseInt(propertyValue);
			dynCorrExcl = (dynCorrExclTime > 0);
			if (dynCorrExcl) {
				mvkeSource.setDynamicCorrelationExclusion(dynCorrExclTime);
				mvkeDest.setDynamicCorrelationExclusion(dynCorrExclTime);
				mvkeJoint.setDynamicCorrelationExclusion(dynCorrExclTime);
			} else {
				mvkeSource.clearDynamicCorrelationExclusion();
				mvkeDest.clearDynamicCorrelationExclusion();
				mvkeJoint.clearDynamicCorrelationExclusion();
			}
		} else if (propertyName.equalsIgnoreCase(FORCE_KERNEL_COMPARE_TO_ALL)) {
			forceCompareToAll = Boolean.parseBoolean(propertyValue);
			mvkeSource.setForceCompareToAll(forceCompareToAll);
			mvkeDest.setForceCompareToAll(forceCompareToAll);
			mvkeJoint.setForceCompareToAll(forceCompareToAll);
		} else {
			// No property was set here
			propertySet = false;
			// try the superclass:
			super.setProperty(propertyName, propertyValue);
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	/**
	 * 
	 * @return the kernel width in use in the calculator
	 */
	public double getKernelWidth() {
		return kernelWidth;
	}
	
	/**
	 * Clone the object - note: while it does create new cloned instances of
	 *  the {@link KernelEstimatorMultiVariate} objects, I think these only
	 *  have shallow copies to the data.
	 * This is enough though to maintain the structure across 
	 *  various {@link #computeSignificance(int)} calls.
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	protected Object clone() throws CloneNotSupportedException {
		MutualInfoCalculatorMultiVariateKernel theClone = 
				(MutualInfoCalculatorMultiVariateKernel) super.clone();
		// Now assign clones of the KernelEstimatorMultiVariate objects:
		theClone.mvkeSource =
				(KernelEstimatorMultiVariate) mvkeSource.clone();
		theClone.mvkeDest =
				(KernelEstimatorMultiVariate) mvkeDest.clone();
		theClone.mvkeJoint =
				(KernelEstimatorMultiVariate) mvkeJoint.clone();
		return theClone;
	}

}
