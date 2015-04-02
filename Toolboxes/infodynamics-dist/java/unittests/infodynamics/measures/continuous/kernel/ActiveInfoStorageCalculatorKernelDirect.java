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

import infodynamics.measures.continuous.ActiveInfoStorageCalculator;
import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;

/**
 * 
 * <p>
 * Implements an active information storage calculator using kernel estimation.
 * </p> 
 * 
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct</li>
 * 		<li>SetProperty() for each property</li>
 *		<li>intialise()</li>
 * 		<li>setObservations(), or [startAddObservations(), addObservations()*, finaliseAddObservations()]
 *   Note: If not using setObservations(), the results from computeLocal or getSignificance
 *    are not likely to be particularly sensible.</li> 
 * 		<li>computeAverageLocalOfObservations() or ComputeLocalOfPreviousObservations()</li>
 * 	</ol>
 * </p>
 * 
 * <p>
 * TODO Use only a single kernel estimator class for the joint space, and compute other
 *  probabilities from this. This will save much time.
 * </p>
 * 
 * <p>
 * Matched against Oliver Obst's c++ implementation on 23/4/2010.
 * </p>
 * 
 * @author Joseph Lizier
 * @see ActiveInfoStorageCalculator
 * @see ActiveInfoStorageCalculatorCorrelationIntegrals
 *
 */
public class ActiveInfoStorageCalculatorKernelDirect
	extends ActiveInfoStorageCalculatorCorrelationIntegrals {

	protected MutualInfoCalculatorMultiVariateKernel miKernel = null;

	// Keep joint vectors so we don't need to regenerate them
	protected double[][] destNextVectors;
	protected double[][] destPastVectors;
	
	/**
	 * Default value for epsilon
	 */
	public static final double DEFAULT_EPSILON = 0.25;
	/**
	 * Kernel width
	 */
	private double epsilon = DEFAULT_EPSILON;

	/**
	 * Creates a new instance of the kernel-estimate style transfer entropy calculator
	 *
	 */
	public ActiveInfoStorageCalculatorKernelDirect() {
		super();
		miKernel = new MutualInfoCalculatorMultiVariateKernel();
	}

	/**
	 * Initialises the calculator with the existing values for k and epsilon
	 * 
	 */
	public void initialise() throws Exception {
		initialise(k, epsilon);
	}

	public void initialise(int k, int tau) throws Exception {
		if (tau != 1) {
			throw new Exception("Not implemented yet for tau != 1");
		}
		initialise(k);
	}
	/**
	 * Initialises the calculator with the existing value for epsilon
	 * 
	 * @param k history length
	 */
	public void initialise(int k) throws Exception {
		this.k = k;
		initialise(k, epsilon);
	}
	
	/**
	 * Initialises the calculator
	 * 
	 * @param k history length
	 * @param epsilon kernel width
	 */
	public void initialise(int k, double epsilon) throws Exception {
		super.initialise(k);
		this.epsilon = epsilon;
		miKernel.initialise(k, 1, epsilon);
		destPastVectors = null;
		destNextVectors = null;
	}

	/**
	 * Set properties for the calculator.
	 * Can include any of the accepted values for
	 * {@link MutualInfoCalculatorMultiVariateKernel#setProperty(String, String)}
	 * except {@link infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF}
	 * 
	 * @param propertyName
	 * @param propertyValue
	 * @throws Exception 
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		super.setProperty(propertyName, propertyValue);
		if (propertyName.equalsIgnoreCase(MutualInfoCalculatorMultiVariateKernel.EPSILON_PROP_NAME)) {
			// Grab epsilon here - we need it for the initialisation:
			epsilon = Double.parseDouble(propertyValue);
			// It will get passed onto miKernel in our initialisation routine
			if (debug) {
				System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
			}
		} if (propertyName.equalsIgnoreCase(MutualInfoCalculatorMultiVariate.PROP_TIME_DIFF)) {
			throw new Exception("Cannot set " + MutualInfoCalculatorMultiVariate.PROP_TIME_DIFF
					+ " property on the ActiveInfoStorageCalculator");
		} else {
			// Pass any other property through to the miKernel object
			miKernel.setProperty(propertyName, propertyValue);
		}
	}

	/**
	 * Flag that the observations are complete, probability distribution functions can now be built.
	 *
	 */
	public void finaliseAddObservations() {
		// First work out the size to allocate the joint vectors, and do the allocation:
		totalObservations = 0;
		for (double[] currentObservation : vectorOfObservations) {
			totalObservations += currentObservation.length - k;
		}
		destPastVectors = new double[totalObservations][k];
		destNextVectors = new double[totalObservations][1];
		
		// Construct the joint vectors from the given observations
		int startObservation = 0;
		for (double[] currentObservation : vectorOfObservations) {
			double[][] currentDestPastVectors = makeJointVectorForPast(currentObservation);
			MatrixUtils.arrayCopy(currentDestPastVectors, 0, 0,
					destPastVectors, startObservation, 0, currentDestPastVectors.length, k);
			double[][] currentDestNextVectors = makeJointVectorForNext(currentObservation);
			MatrixUtils.arrayCopy(currentDestNextVectors, 0, 0,
					destNextVectors, startObservation, 0, currentDestNextVectors.length, 1);
			startObservation += currentObservation.length - k;
		}
		
		// Now set the joint vectors in the kernel estimator
		try {
			miKernel.setObservations(destPastVectors, destNextVectors);
		} catch (Exception e) {
			// Should only throw where our vector sizes don't match
			throw new RuntimeException(e);
		}

		// Store whether there was more than one observation set:
		addedMoreThanOneObservationSet = vectorOfObservations.size() > 1;
		
		// And clear the vector of observations
		vectorOfObservations = null;
	}
	
	/**
	 * <p>Computes the average Active Info Storage for the previously supplied observations</p> 
	 * 
	 */
	public double computeAverageLocalOfObservations() throws Exception {
		lastAverage = miKernel.computeAverageLocalOfObservations();
		return lastAverage;
	}

	/**
	 * Computes the local active information storage for the previous supplied observations.
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
		double[] local = miKernel.computeLocalOfPreviousObservations();
		lastAverage = miKernel.getLastAverage();
		if (!addedMoreThanOneObservationSet) {
			double[] localsToReturn = new double[local.length + k];
			System.arraycopy(local, 0, localsToReturn, k, local.length);
			return localsToReturn;
		} else {
			return local;
		}
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
		return miKernel.computeSignificance(numPermutationsToCheck);
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
		
		return miKernel.computeSignificance(newOrderings);
	}

	public void setDebug(boolean debug) {
		super.setDebug(debug);
		miKernel.setDebug(debug);
	}

	public double[] computeLocalUsingPreviousObservations(
			double[] newObservations) throws Exception {
		throw new Exception("Not yet implemented");
	}
}
