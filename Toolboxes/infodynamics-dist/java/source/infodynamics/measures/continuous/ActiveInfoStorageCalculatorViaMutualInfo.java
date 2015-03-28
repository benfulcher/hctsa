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

package infodynamics.measures.continuous;

import java.util.Vector;

import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;

/**
 * An Active Information Storage (AIS) calculator (implementing {@link ActiveInfoStorageCalculator})
 * which is affected using a 
 * given Mutual Information (MI) calculator (implementing
 * {@link MutualInfoCalculatorMultiVariate}) to make the calculations.
 * 
 * <p>Usage is as per the paradigm outlined for {@link ActiveInfoStorageCalculator},
 * except that in the constructor(s) for this class the implementation for
 * a {@link MutualInfoCalculatorMultiVariate} must be supplied.
 * </p>
 * 
 * <p>This class <i>may</i> be used directly, however users are advised that
 * several child classes are available which already plug-in the various MI estimators
 * to provide AIS calculators (taking specific caution associated with
 * each type of estimator):</p>
 * <ul>
 * 	<li>{@link infodynamics.measures.continuous.gaussian.ActiveInfoStorageCalculatorGaussian}</li>
 * 	<li>{@link infodynamics.measures.continuous.kernel.ActiveInfoStorageCalculatorKernel}</li>
 * 	<li>{@link infodynamics.measures.continuous.kraskov.ActiveInfoStorageCalculatorKraskov}</li>
 * </ul>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>J.T. Lizier, M. Prokopenko and A.Y. Zomaya,
 * 		<a href="http://dx.doi.org/10.1016/j.ins.2012.04.016">
 * 		"Local measures of information storage in complex distributed computation"</a>,
 * 		Information Sciences, vol. 208, pp. 39-54, 2012.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @see ActiveInfoStorageCalculator
 */
public class ActiveInfoStorageCalculatorViaMutualInfo implements
		ActiveInfoStorageCalculator {

	/**
	 * The underlying mutual information calculator
	 */
	protected MutualInfoCalculatorMultiVariate miCalc;
	/**
	 * Length of past history to consider (embedding length)
	 */
	protected int k = 1;
	/**
	 * Embedding delay to use between elements of the embeding vector.
	 * We're hard-coding a delay of 1 between the history vector and the next 
	 *  observation however.
	 */
	protected int tau = 1;
	/**
	 * Whether debug mode is on
	 */
	protected boolean debug = false;
	
	/**
	 * Construct using an instantiation of the named MI calculator
	 * 
	 * @param miCalculatorClassName fully qualified class name of the MI calculator to instantiate
	 * @throws InstantiationException
	 * @throws IllegalAccessException
	 * @throws ClassNotFoundException
	 */
	public ActiveInfoStorageCalculatorViaMutualInfo(String miCalculatorClassName)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		@SuppressWarnings("unchecked")
		Class<MutualInfoCalculatorMultiVariate> miClass = 
				(Class<MutualInfoCalculatorMultiVariate>) Class.forName(miCalculatorClassName);
		MutualInfoCalculatorMultiVariate miCalc = miClass.newInstance();
		construct(miCalc);
	}
	
	/**
	 * Construct using an instantiation of the given MI class
	 * 
	 * @param miCalcClass Class of the MI calculator to instantiate and use
	 * @throws InstantiationException
	 * @throws IllegalAccessException
	 */
	protected ActiveInfoStorageCalculatorViaMutualInfo(Class<MutualInfoCalculatorMultiVariate> miCalcClass)
			throws InstantiationException, IllegalAccessException {
		MutualInfoCalculatorMultiVariate miCalc = miCalcClass.newInstance();
		construct(miCalc);
	}
	
	/**
	 * Construct using the given (constructed but not initialised)
	 * MI calculator.
	 * 
	 * @param miCalc MI calculator which is already constructed but
	 *  there has not been a call to its {@link MutualInfoCalculatorMultiVariate#initialise()}
	 *  method yet 
	 */
	protected ActiveInfoStorageCalculatorViaMutualInfo(MutualInfoCalculatorMultiVariate miCalc) {
		construct(miCalc);
	}

	/**
	 * Internal routine to execute common code for constructing an instance
	 * 
	 * @param miCalc
	 */
	protected void construct(MutualInfoCalculatorMultiVariate miCalc) {
		this.miCalc = miCalc;
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#initialise()
	 */
	@Override
	public void initialise() throws Exception {
		initialise(k, tau); // Initialise with current value of k
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#initialise(int)
	 */
	@Override
	public void initialise(int k) throws Exception {
		initialise(k, tau);
	}

	/**
	 * {@inheritDoc}
	 * 
	 * <p>All child classes <b>must</b> call this routine on this as the super class
	 *  once they have finished executing their specialised code
	 *  for their {@link #initialise()} implementations.
	 * </p>
	 * 
	 */
	@Override
	public void initialise(int k, int tau) throws Exception {
		this.k = k;
		this.tau = tau;
		miCalc.initialise(k, 1);
	}

	/**
	 * Set properties for the underlying calculator implementation.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 *  
	 * <p>Allowable property names include:</p>
	 * <ul>
	 * 	<li>Those defined for the {@link ActiveInfoStorageCalculator} interface
	 *      (i.e. {@link #K_PROP_NAME} or {@link #TAU_PROP_NAME})</li>
	 *  <li>Any properties defined for the underlying
	 *     {@link MutualInfoCalculatorMultiVariate#setProperty(String, String)} implementation,
	 *     <b>however</b> the user is <b>not</b> allowed to set the property 
	 *     {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF} here.
	 *     This would set a time difference from the history vector to the next
	 *     step, which we currently do not allow.
	 *     (If we change our mind one day and allow it, we could implement
	 *     it simply by letting the time diff property be set here).</li>
	 * </ul>
	 *  
	 * <p>Note that implementing classes may defined additional properties.</p>
	 *
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		if (propertyName.equalsIgnoreCase(MutualInfoCalculatorMultiVariate.PROP_TIME_DIFF)) {
			throw new Exception("Cannot set " + MutualInfoCalculatorMultiVariate.PROP_TIME_DIFF
					+ " property on the ActiveInfoStorageCalculator");
		}
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(K_PROP_NAME)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(TAU_PROP_NAME)) {
			tau = Integer.parseInt(propertyValue);
		} else {
			// No property was set on this class, assume it is for the underlying
			//  MI calculator
			miCalc.setProperty(propertyName, propertyValue);
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#setObservations(double[])
	 */
	@Override
	public void setObservations(double[] observations) throws Exception {
		if (observations.length - (k-1)*tau - 1 <= 0) {
			// There are no observations to add here
			throw new Exception("Not enough observations to set here given k and tau");
		}
		double[][] currentDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(observations, k, tau, (k-1)*tau, observations.length - (k-1)*tau - 1);
		double[][] currentDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(observations, 1, (k-1)*tau + 1, observations.length - (k-1)*tau - 1);
		miCalc.setObservations(currentDestPastVectors, currentDestNextVectors);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#startAddObservations()
	 */
	@Override
	public void startAddObservations() {
		miCalc.startAddObservations();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#addObservations(double[])
	 */
	@Override
	public void addObservations(double[] observations) throws Exception {
		if (observations.length - (k-1)*tau - 1 <= 0) {
			// There are no observations to add here
			// Don't throw an exception, do nothing since more observations
			//  can be added later.
			return;
		}
		double[][] currentDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(observations, k, tau, (k-1)*tau, observations.length - (k-1)*tau - 1);
		double[][] currentDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(observations, 1, (k-1)*tau + 1, observations.length - (k-1)*tau - 1);
		miCalc.addObservations(currentDestPastVectors, currentDestNextVectors);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#addObservations(double[], int, int)
	 */
	@Override
	public void addObservations(double[] observations, int startTime,
			int numTimeSteps) throws Exception {
		addObservations(MatrixUtils.select(observations, startTime, numTimeSteps));
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#finaliseAddObservations()
	 */
	@Override
	public void finaliseAddObservations() throws Exception {
		miCalc.finaliseAddObservations();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#setObservations(double[], boolean[])
	 */
	@Override
	public void setObservations(double[] observations, boolean[] valid)
			throws Exception {
		Vector<int[]> startAndEndTimePairs = computeStartAndEndTimePairs(valid);
		
		// We've found the set of start and end times for this pair
		startAddObservations();
		for (int[] timePair : startAndEndTimePairs) {
			int startTime = timePair[0];
			int endTime = timePair[1];
			addObservations(observations, startTime, endTime - startTime + 1);
		}
		finaliseAddObservations();
	}

	/**
	 * Compute a vector of start and end pairs of time points, between which we have
	 *  valid series of observations.
	 * 
	 * <p>This method is made public so it can be used if one wants to compute the number of
	 *  observations prior to making a call to {@link #setObservations(double[], boolean[])}.</p>
	 * 
	 * @param valid a time series (with indices the same as observations)
	 *  indicating whether the entry in observations at that index is valid; 
	 *  we only take vectors as samples to add to the observation set where
	 *  all points in the time series (even between points in 
	 *  the embedded k-vector with embedding delays) are valid.
	 * @return a vector for start and end time pairs of valid series
	 *  of observations (as defined by <code>valid</code>).
	 */
	public Vector<int[]> computeStartAndEndTimePairs(boolean[] valid) {
		// Scan along the data avoiding invalid values
		int startTime = 0;
		int endTime = 0;
		boolean lookingForStart = true;
		Vector<int[]> startAndEndTimePairs = new Vector<int[]>();
		for (int t = 0; t < valid.length; t++) {
			if (lookingForStart) {
				// Precondition: startTime holds a candidate start time
				if (valid[t]) {
					// This point is OK at the destination
					if (t - startTime < (k-1)*tau+1) {
						// We're still checking the past history only, so
						continue;
					} else {
						// We've got the full past history ok
						// set a candidate endTime
						endTime = t;
						lookingForStart = false;
						if (t == valid.length - 1) {
							// we need to terminate now
							int[] timePair = new int[2];
							timePair[0] = startTime;
							timePair[1] = endTime;
							startAndEndTimePairs.add(timePair);
							// System.out.printf("t_s=%d, t_e=%d\n", startTime, endTime);
						}
					}
				} else {
					// We need to keep looking.
					// Move the potential start time to the next point
					startTime = t + 1;
				}
			} else {
				// Precondition: startTime holds the start time for this set, 
				//  endTime holds a candidate end time
				// Check if we can include the current time step
				boolean terminateSequence = false;
				if (valid[t]) {
					// We can extend
					endTime = t;
				} else {
					terminateSequence = true;
				}
				if (t == valid.length - 1) {
					// we need to terminate the sequence anyway
					terminateSequence = true;
				}
				if (terminateSequence) {
					// This section is done
					int[] timePair = new int[2];
					timePair[0] = startTime;
					timePair[1] = endTime;
					startAndEndTimePairs.add(timePair);
					// System.out.printf("t_s=%d, t_e=%d\n", startTime, endTime);
					lookingForStart = true;
					startTime = t + 1;
				}
			}
		}
		return startAndEndTimePairs;
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#computeAverageLocalOfObservations()
	 */
	@Override
	public double computeAverageLocalOfObservations() throws Exception {
		return miCalc.computeAverageLocalOfObservations();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#computeLocalOfPreviousObservations()
	 */
	@Override
	public double[] computeLocalOfPreviousObservations() throws Exception {
		double[] local = miCalc.computeLocalOfPreviousObservations();
		if (!miCalc.getAddedMoreThanOneObservationSet()) {
			double[] localsToReturn = new double[local.length + (k-1)*tau + 1];
			System.arraycopy(local, 0, localsToReturn, (k-1)*tau + 1, local.length);
			return localsToReturn;
		} else {
			return local;
		}
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#computeLocalUsingPreviousObservations(double[])
	 */
	@Override
	public double[] computeLocalUsingPreviousObservations(double[] newObservations) throws Exception {
		if (newObservations.length - (k-1)*tau - 1 <= 0) {
			// There are no observations to compute for here
			return new double[newObservations.length];
		}
		double[][] newDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(newObservations, k, tau, (k-1)*tau, newObservations.length - (k-1)*tau - 1);
		double[][] newDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(newObservations, 1, (k-1)*tau + 1, newObservations.length - (k-1)*tau - 1);
		double[] local = miCalc.computeLocalUsingPreviousObservations(newDestPastVectors, newDestNextVectors);
		// Pad the front of the array with zeros where local AIS isn't defined:
		double[] localsToReturn = new double[local.length + (k-1)*tau + 1];
		System.arraycopy(local, 0, localsToReturn, (k-1)*tau + 1, local.length);
		return localsToReturn;

	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#computeSignificance(int)
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		return miCalc.computeSignificance(numPermutationsToCheck);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#computeSignificance(int[][])
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception {
		return miCalc.computeSignificance(newOrderings);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#setDebug(boolean)
	 */
	@Override
	public void setDebug(boolean debug) {
		this.debug = debug;
		miCalc.setDebug(debug);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#getLastAverage()
	 */
	@Override
	public double getLastAverage() {
		return miCalc.getLastAverage();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ActiveInfoStorageCalculator#getNumObservations()
	 */
	@Override
	public int getNumObservations() throws Exception {
		return miCalc.getNumObservations();
	}
}
