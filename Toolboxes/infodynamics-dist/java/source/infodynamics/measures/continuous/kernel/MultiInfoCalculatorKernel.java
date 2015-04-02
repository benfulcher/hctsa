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

import infodynamics.measures.continuous.MultiInfoCalculator;
import infodynamics.measures.continuous.MultiInfoCalculatorCommon;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential multi-information of a given multivariate set of
 *  observations (implementing {@link MultiInfoCalculator}),
 *  using box-kernel estimation.
 *  For details on box-kernel estimation, see Kantz and Schreiber (below).</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link MultiInfoCalculator},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #MultiInfoCalculatorKernel()}.</li>
 *  <li>Further properties are available, see {@link #setProperty(String, String)};</li>
 *  <li>An additional {@link #initialise(int, double)} option;</li>
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
public class MultiInfoCalculatorKernel 
	extends MultiInfoCalculatorCommon {

	/**
	 * Marginal kernel density PDF estimators
	 */
	protected KernelEstimatorUniVariate[] svkeMarginals = null;
	/**
	 * Joint space kernel density PDF estimator
	 */
	protected KernelEstimatorMultiVariate mvkeJoint = null;


	private boolean dynCorrExcl = false;
	private int dynCorrExclTime = 100;
	/**
	 * Property name for a dynamics exclusion time window (see Kantz and Schreiber),
	 * default is 0 which means no dynamic exclusion window.
	 */
	public static final String DYN_CORR_EXCL_TIME_NAME = "DYN_CORR_EXCL";

	/**
	 * Property name for the kernel width
	 */
	public static final String KERNEL_WIDTH_PROP_NAME = "KERNEL_WIDTH";
	/**
	 * Legacy property name for the kernel width
	 */
	public static final String EPSILON_PROP_NAME = "EPSILON";
	/**
	 * Default value for kernel width
	 */
	public static final double DEFAULT_KERNEL_WIDTH = 0.25;
	/**
	 * Kernel width currently in use
	 */
	private double kernelWidth = DEFAULT_KERNEL_WIDTH;

	/**
	 * Construct an instance
	 */
	public MultiInfoCalculatorKernel() {
		mvkeJoint = new KernelEstimatorMultiVariate();
	}

	@Override
	public void initialise(int dimensions) {
		// Super class' initialise() will be called
		//  from the following:
		initialise(dimensions, kernelWidth);
	}

	/**
	 * Initialise the calculator for (re-)use, with a specific kernel width and
	 * with number of joint variables specified, and existing
	 * (or default) values of other parameters,.
	 * Clears an PDFs of previously supplied observations.
	 *
	 * @param dimensions the number of joint variables
	 * @param kernelWidth if {@link #PROP_NORMALISE} property has
	 *  been set, then this kernel width corresponds to the number of
	 *  standard deviations from the mean (otherwise it is an absolute value)
	 */
	public void initialise(int dimensions, double epsilon) {
		this.kernelWidth = epsilon;
		if (this.dimensions != dimensions) {
			// Need to create a new array of marginal kernel estimators 
			this.dimensions = dimensions;
			svkeMarginals = new KernelEstimatorUniVariate[dimensions];
			for (int i = 0; i < dimensions; i++) {
				svkeMarginals[i] = new KernelEstimatorUniVariate();
				svkeMarginals[i].setNormalise(normalise);
				if (dynCorrExcl) {
					svkeMarginals[i].setDynamicCorrelationExclusion(dynCorrExclTime);
				} else {
					svkeMarginals[i].clearDynamicCorrelationExclusion();
				}
			}
		}
		// Initialise the marginal kernel estimators
		for (int i = 0; i < dimensions; i++) {
			svkeMarginals[i].initialise(epsilon);
		}
		// Initialise the joint kernel estimator
		mvkeJoint.initialise(dimensions, epsilon);
		// Now call the super class to handle the common variables:
		super.initialise(dimensions);
	}

	/**
	 * <p>Set properties for the kernel multi-information calculator.
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
	 * 		<li>{@link #DYN_CORR_EXCL_TIME_NAME} -- a dynamics exclusion time window (see Kantz and Schreiber),
	 * 			default is 0 which means no dynamic exclusion window.</li>
	 *  	<li>any valid properties for {@link MultiInfoCalculatorCommon#setProperty(String, String)}.</li>
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
		} else if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
			super.setProperty(propertyName, propertyValue);
			// More to do with this property locally:
			for (int d = 0; d < dimensions; d++) {
				svkeMarginals[d].setNormalise(normalise);
			}
			mvkeJoint.setNormalise(normalise);
		} else if (propertyName.equalsIgnoreCase(DYN_CORR_EXCL_TIME_NAME)) {
			dynCorrExclTime = Integer.parseInt(propertyValue);
			dynCorrExcl = (dynCorrExclTime > 0);
			if (dynCorrExcl) {
				for (int d = 0; d < dimensions; d++) {
					svkeMarginals[d].setDynamicCorrelationExclusion(dynCorrExclTime);
				}
				mvkeJoint.setDynamicCorrelationExclusion(dynCorrExclTime);
			} else {
				for (int d = 0; d < dimensions; d++) {
					svkeMarginals[d].clearDynamicCorrelationExclusion();
				}
				mvkeJoint.clearDynamicCorrelationExclusion();
			}
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

	@Override
	public void startAddObservations() {
		if (dynCorrExcl) {
			// We have not properly implemented dynamic correlation exclusion for
			//  multiple observation sets, so throw an error
			throw new RuntimeException("Addition of multiple observation sets is not currently " +
					"supported with property DYN_CORR_EXCL set");
		}
		super.startAddObservations();
	}
	
	@Override
	public void finaliseAddObservations() throws Exception {
		super.finaliseAddObservations();

		for (int d = 0; d < dimensions; d++) {
			svkeMarginals[d].setObservations(MatrixUtils.selectColumn(observations, d));
		}
		mvkeJoint.setObservations(observations);
	}

	@Override
	public double computeAverageLocalOfObservations() {
		double mi = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			if (debug) {
				System.out.print(b + ": ");
			}
			double marginalProbProducts = 1.0;
			for (int d = 0; d < dimensions; d++) {
				double marginalProb = svkeMarginals[d].getProbability(observations[b][d], b);
				if (debug) {
					System.out.print(observations[b][d] + " p=" + marginalProb + ", ");
				}
				marginalProbProducts *= marginalProb;
			}
			double probJoint = mvkeJoint.getCount(observations[b], b);
			double logTerm = 0.0;
			double cont = 0.0;
			if (probJoint > 0.0) {
				// TODO Should probably check that marginalProbProducts has not 
				//  gone to zero (this is possible with several multiplications
				//  of 1/N, though is unlikely). I'm not sure what we would do if
				//  it had gone to zero though ... ignore this value?
				logTerm = probJoint / marginalProbProducts;
				cont = Math.log(logTerm);
			}
			mi += cont;
			if (debug) {
				System.out.println(", p(joint) = " + probJoint
						+ " -> " + logTerm + " -> " + (cont/Math.log(2.0)) + " -> sum: " + (mi/Math.log(2.0)));
			}
		}
		lastAverage = mi / (double) totalObservations / Math.log(2.0);
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
			double prob = mvkeJoint.getCount(observations[b], b);
			double cont = 0.0;
			if (prob > 0.0) {
				cont = - Math.log(prob);
			}
			entropy += cont;
			if (debug) {
				System.out.println(b + ": " + prob
						+ " -> " + cont/Math.log(2.0) + " -> sum: " + (entropy/Math.log(2.0)));
			}
		}
		return entropy / (double) totalObservations / Math.log(2.0);
	}

	/**
	 * Extra utility method to return the entropy of the given individual variable
	 * 
	 * @param variableIndex which variable to compute the entropy for
	 * @return entropy of given variable in bits.
	 */
	public double computeAverageMarginalEntropy(int variableIndex) {
		double entropy = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob = svkeMarginals[variableIndex].getProbability(observations[b][variableIndex], b);
			double cont = 0.0;
			if (prob > 0.0) {
				cont = -Math.log(prob);
			}
			entropy += cont;
			if (debug) {
				System.out.println(b + ": " + prob
						+ " -> " + cont/Math.log(2.0) + " -> sum: " + (entropy/Math.log(2.0)));
			}
		}
		return entropy / (double) totalObservations / Math.log(2.0);
	}

	/**
	 * Extra utility method to return the information distance between the
	 * marginal variables (we've generalised the definition from pair-wise here).
	 * 
	 * @return information distance between the source and destination in bits.
	 */
	public double computeAverageInfoDistanceOfObservations() {
		double infoDistance = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double marginalProbProducts = 1.0;
			for (int d = 0; d < dimensions; d++) {
				marginalProbProducts *= svkeMarginals[d].getProbability(observations[b][d], b);
			}
			double probJoint = mvkeJoint.getProbability(observations[b], b);
			double logTerm = 0.0;
			double cont = 0.0;
			if (probJoint > 0.0) {
				// TODO Should probably check that marginalProbProducts has not 
				//  gone to zero (this is possible with several multiplications
				//  of 1/N, though is unlikely). I'm not sure what we would do if
				//  it had gone to zero though ... ignore this value?
				//  It's not easy to ignore since we can't say p log p -> 0 here
				//  because marginalProbProducts is not multiplying out the front
				logTerm = marginalProbProducts / (probJoint * probJoint);
				cont = Math.log(logTerm);
			}
			infoDistance += cont;
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + (cont/Math.log(2.0)) +
						" -> sum: " + (infoDistance/Math.log(2.0)));
			}
		}
		return infoDistance / (double) totalObservations / Math.log(2.0);
	}

	@Override
	public double[] computeLocalOfPreviousObservations() throws Exception {
		return computeLocalUsingPreviousObservations(observations, true);
	}

	/**
	 * Compute the local multi-information values for each of the
	 * supplied samples in <code>states</code>.
	 * 
	 * <p>PDFs are computed using all of the previously supplied
	 * observations, but not those in <code>states</code>
	 * (unless they were
	 * some of the previously supplied samples).</p>
	 * 
	 * <p>Note that calls to this method will not harness
	 *  dynamic correlation exclusion (if set)
	 *  since we don't know whether it's the same time set or not.</p>
	 */
	@Override
	public double[] computeLocalUsingPreviousObservations(double states[][]) {
		return computeLocalUsingPreviousObservations(states, false);
	}
	
	/**
	 * Internal implementation for {@link #computeLocalUsingPreviousObservations(double[][])}
	 *  and {@link #computeLocalOfPreviousObservations()}.
	 *  
	 * @param states series of multivariate observations
	 *  (first index is time or observation index, second is variable number)
	 * @param isOurPreviousObservations true if we are implementing 
	 *  {@link #computeLocalOfPreviousObservations()}, false for 
	 *  {@link #computeLocalUsingPreviousObservations(double[][])}
	 * @return the series of local multi-information values.
	 */
	protected double[] computeLocalUsingPreviousObservations(double states[][],
			boolean isOurPreviousObservations) {
		double mi = 0.0;
		int timeSteps = states.length;
		double[] localMi = new double[timeSteps];
		double probJoint;
		for (int b = 0; b < timeSteps; b++) {
			double marginalProbProducts = 1.0;
			for (int d = 0; d < dimensions; d++) {
				if (isOurPreviousObservations) {
					marginalProbProducts *= svkeMarginals[d].getProbability(states[b][d], b);
				} else {
					marginalProbProducts *= svkeMarginals[d].getProbability(states[b][d]);
				}
			}
			if (isOurPreviousObservations) {
				probJoint = mvkeJoint.getProbability(states[b], b);
			} else {
				probJoint = mvkeJoint.getProbability(states[b]);
			}
			double logTerm = 0.0;
			localMi[b] = 0.0;
			if (probJoint > 0.0) {
				// TODO Should probably check that marginalProbProducts has not 
				//  gone to zero (this is possible with several multiplications
				//  of 1/N, though is unlikely). I'm not sure what we would do if
				//  it had gone to zero though ... ignore this value?
				logTerm = probJoint / marginalProbProducts;
				localMi[b] = Math.log(logTerm) / Math.log(2.0);
			}
			mi += localMi[b];
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + localMi[b] + " -> sum: " + mi);
			}
		}
		lastAverage = mi / (double) totalObservations;
		return localMi;
	}

	
	/**
	 * Compute the local joint entropy values of the previously provided
	 *  observations.
	 *  
	 * @return the local joint entropies in bits
	 */
	public double[] computeLocalJointEntropyOfPreviousObservations() throws Exception {
		return computeLocalJointEntropyUsingPreviousObservations(observations, true);
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
	public double[] computeLocalJointEntropyUsingPreviousObservations(double states[][]) {
		return computeLocalJointEntropyUsingPreviousObservations(states, false);
	}
	
	/**
	 * Internal implementation for {@link #computeLocalJointEntropyUsingPreviousObservations(double[][])}
	 *  and {@link #computeLocalJointEntropyOfPreviousObservations()}.
	 *  
	 * @param states series of multivariate observations
	 *  (first index is time or observation index, second is variable number)
	 * @param isOurPreviousObservations true if we are implementing 
	 *  {@link #computeLocalJointEntropyOfPreviousObservations()}, false for 
	 *  {@link #computeLocalJointEntropyUsingPreviousObservations(double[][])}
	 * @return the series of local joint entropy values.
	 */
	protected double[] computeLocalJointEntropyUsingPreviousObservations(
			double states[][], boolean isOurPreviousObservations) {

		int timeSteps = states.length;
		double[] localJoint = new double[timeSteps];
		double probJoint;
		for (int b = 0; b < totalObservations; b++) {
			if (isOurPreviousObservations) {
				probJoint = mvkeJoint.getProbability(states[b], b);
			} else {
				probJoint = mvkeJoint.getProbability(states[b]);
			}
			localJoint[b] = 0.0;
			if (probJoint > 0.0) {
				localJoint[b] = - Math.log(probJoint) / Math.log(2.0);
			}
			if (debug) {
				System.out.println(b + ": " + probJoint + " -> " + localJoint[b]);
			}
		}
		return localJoint;
	}

	/**
	 * Compute the local entropy values for the previously provided
	 *  observations for the given variable 
	 *  (using those previous observations to compute the PDFs).
	 * 
	 * @param variableIndex which variable to compute local entropies for
	 * @return array of local entropies for the given variable
	 */
	public double[] computeLocalMarginalEntropyOfPreviousObservations(int variableIndex) {
		return computeLocalMarginalEntropyUsingPreviousObservations(observations, variableIndex, true);
	}

	/**
	 * Compute the local entropy values for the given variable, for the 
	 * give observations,
	 * using the previously provided
	 *  observations to compute the probabilities.
	 * Calls to this method will not harness dynamic correlation exclusion (if set)
	 *  since we don't know whether it's the same time set or not.
	 * 
	 * @param states provided observations
	 * @param variableIndex which variable to compute local entropies for
	 * @return array of local entropies for these observations for this variable
	 */
	public double[] computeLocalMarginalEntropyUsingPreviousObservations(double states[][], int variableIndex) {
		return computeLocalMarginalEntropyUsingPreviousObservations(states, variableIndex, false);
	}
	
	/**
	 * Internal implementation for
	 * {@link #computeLocalMarginalEntropyUsingPreviousObservations(double[][], int)}
	 *  and {@link #computeLocalMarginalEntropyOfPreviousObservations(int)}.
	 *  
	 * @param states series of multivariate observations
	 *  (first index is time or observation index, second is variable number)
	 * @param variableIndex which variable to compute local entropies for
	 * @param isOurPreviousObservations true if we are implementing 
	 *  {@link #computeLocalMarginalEntropyOfPreviousObservations(int)}, false for 
	 *  {@link #computeLocalMarginalEntropyUsingPreviousObservations(double[][], int)}
	 * @return the series of local marginal entropy values for the given variable.
	 */
	protected double[] computeLocalMarginalEntropyUsingPreviousObservations(
			double states[][], int variableIndex, boolean isOurPreviousObservations) {
		int timeSteps = states.length;
		double[] localEntropy = new double[timeSteps];
		double prob;
		for (int b = 0; b < totalObservations; b++) {
			if (isOurPreviousObservations) {
				prob = svkeMarginals[variableIndex].getProbability(states[b][variableIndex], b);
			} else {
				prob = svkeMarginals[variableIndex].getProbability(states[b][variableIndex]);
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
		return computeLocalInfoDistanceUsingPreviousObservations(observations, true);
	}

	/**
	 * Compute the local Info distance values for these given values, using the previously provided
	 *  observations to compute the probabilities.
	 * Calls to this method will not harness dynamic correlation exclusion (if set)
	 *  since we don't know whether it's the same time set or not.
	 * 
	 * @return array of local information distances.
	 */
	public double[] computeLocalInfoDistanceUsingPreviousObservations(double[][] states) {
		return computeLocalInfoDistanceUsingPreviousObservations(states, false);
	}

	/**
	 * Internal implementation for
	 * {@link #computeLocalInfoDistanceUsingPreviousObservations(double[][])}
	 *  and {@link #computeLocalInfoDistanceOfPreviousObservations()}.
	 *  
	 * @param states series of multivariate observations
	 *  (first index is time or observation index, second is variable number)
	 * @param isOurPreviousObservations true if we are implementing 
	 *  {@link #computeLocalInfoDistanceOfPreviousObservations()}, false for 
	 *  {@link #computeLocalInfoDistanceUsingPreviousObservations(double[][])}
	 * @return the series of local info distance values.
	 */
	protected double[] computeLocalInfoDistanceUsingPreviousObservations(
			double[][] states, boolean isOurPreviousObservations) {
		
		int timeSteps = states.length;
		double[] localInfoDistance = new double[timeSteps];
		double probJoint;
		for (int b = 0; b < timeSteps; b++) {
			double marginalProbProducts = 1.0;
			for (int d = 0; d < dimensions; d++) {
				if (isOurPreviousObservations) {
					marginalProbProducts *= svkeMarginals[d].getProbability(states[b][d], b);
				} else {
					marginalProbProducts *= svkeMarginals[d].getProbability(states[b][d]);
				}
			}
			if (isOurPreviousObservations) {
				probJoint = mvkeJoint.getProbability(states[b], b);
			} else {
				probJoint = mvkeJoint.getProbability(states[b]);
			}
			double logTerm = 0.0;
			localInfoDistance[b] = 0.0;
			if (probJoint > 0.0) {
				// TODO Should probably check that marginalProbProducts has not 
				//  gone to zero (this is possible with several multiplications
				//  of 1/N, though is unlikely). I'm not sure what we would do if
				//  it had gone to zero though ... ignore this value?
				//  It's not easy to ignore since we can't say p log p -> 0 here
				//  because marginalProbProducts is not multiplying out the front
				logTerm = marginalProbProducts / (probJoint * probJoint);
				localInfoDistance[b] = Math.log(logTerm) / Math.log(2.0);
			}
			if (debug) {
				System.out.println(b + ": " + logTerm + " -> " + localInfoDistance[b]);
			}
		}
		return localInfoDistance;
	}
}
