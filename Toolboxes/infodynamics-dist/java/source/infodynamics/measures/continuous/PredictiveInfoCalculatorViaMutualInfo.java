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
 * A Predictive Information (PI) / Excess Entropy calculator
 * (implementing {@link PredictiveInfoCalculator})
 * which is affected using a 
 * given Mutual Information (MI) calculator (implementing
 * {@link MutualInfoCalculatorMultiVariate}) to make the calculations.
 * 
 * <p>Usage is as per the paradigm outlined for {@link PredictiveInfoCalculator},
 * except that in the constructor(s) for this class the implementation for
 * a {@link MutualInfoCalculatorMultiVariate} must be supplied.
 * </p>
 * 
 * <p>This class <i>may</i> be used directly, however users are advised that
 * several child classes are available which already plug-in the various MI estimators
 * to provide AIS calculators (taking specific caution associated with
 * each type of estimator):</p>
 * <ul>
 * 	<li>{@link infodynamics.measures.continuous.gaussian.PredictiveInfoCalculatorGaussian}</li>
 * 	<li>{@link infodynamics.measures.continuous.kernel.PredictiveInfoCalculatorKernel}</li>
 * 	<li>{@link infodynamics.measures.continuous.kraskov.PredictiveInfoCalculatorKraskov}</li>
 * </ul>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Bialek, W., Nemenman, I., and Tishby, N.,
 *  <a href="http://dx.doi.org/10.1016/S0378-4371(01)00444-7">
 * 	"Complexity through nonextensivity"</a>,
 *  Physica A, 302, 89-99. (2001).</li>
 * 	<li>J. P. Crutchfield, D. P. Feldman,
 *  <a href="http://dx.doi.org/10.1063/1.1530990">
 * 	"Regularities Unseen, Randomness Observed: Levels of Entropy Convergence"</a>,
 *  Chaos, Vol. 13, No. 1. (2003), pp. 25-54.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @see ActiveInfoStorageCalculator
 */
public class PredictiveInfoCalculatorViaMutualInfo implements
		PredictiveInfoCalculator {

	/**
	 * The underlying mutual information calculator
	 */
	protected MutualInfoCalculatorMultiVariate miCalc;
	/**
	 * Length of past and future vectors to consider (embedding length)
	 */
	protected int k = 1;
	/**
	 * Embedding delay to use between elements of the embedding vector.
	 * We're hard-coding a delay of 1 between the history vector and the future
	 *  vector however.
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
	public PredictiveInfoCalculatorViaMutualInfo(String miCalculatorClassName)
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
	protected PredictiveInfoCalculatorViaMutualInfo(Class<MutualInfoCalculatorMultiVariate> miCalcClass)
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
	protected PredictiveInfoCalculatorViaMutualInfo(MutualInfoCalculatorMultiVariate miCalc) {
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
	 * @see infodynamics.measures.continuous.PredictiveInfoCalculator#initialise()
	 */
	@Override
	public void initialise() throws Exception {
		initialise(k, tau); // Initialise with current value of k
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.PredictiveInfoCalculator#initialise(int)
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
		miCalc.initialise(k, k);
	}

	/**
	 * Set properties for the underlying calculator implementation.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 *  
	 * <p>Allowable property names include:</p>
	 * <ul>
	 * 	<li>Those defined for the {@link PredictiveInfoCalculator} interface
	 *      (i.e. {@link #K_PROP_NAME} or {@link K_EMBEDDING_PROP_NAME}
	 *      or {@link #TAU_PROP_NAME})</li>
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
					+ " property on the PredictiveInfoCalculator");
		}
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(K_PROP_NAME) ||
			propertyName.equalsIgnoreCase(K_EMBEDDING_PROP_NAME)) {
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
	 * @see infodynamics.measures.continuous.PredictiveInfoCalculator#setObservations(double[])
	 */
	@Override
	public void setObservations(double[] observations) throws Exception {
		if (observations.length - 2*(k-1)*tau - 1 <= 0) {
			// There are no observations to add here
			throw new Exception("Not enough observations to set here given k and tau");
		}
		double[][] currentPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(observations, k, tau, (k-1)*tau, observations.length - 2*(k-1)*tau - 1);
		double[][] currentNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(observations, k, 2*(k-1)*tau + 1, observations.length - 2*(k-1)*tau - 1);
		miCalc.setObservations(currentPastVectors, currentNextVectors);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.PredictiveInfoCalculator#startAddObservations()
	 */
	@Override
	public void startAddObservations() {
		miCalc.startAddObservations();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.PredictiveInfoCalculator#addObservations(double[])
	 */
	@Override
	public void addObservations(double[] observations) throws Exception {
		if (observations.length - 2*(k-1)*tau - 1 <= 0) {
			// There are no observations to add here
			// Don't throw an exception, do nothing since more observations
			//  can be added later.
			return;
		}
		double[][] currentPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(observations, k, tau, (k-1)*tau, observations.length - 2*(k-1)*tau - 1);
		double[][] currentNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(observations, 1, 2*(k-1)*tau + 1, observations.length - 2*(k-1)*tau - 1);
		miCalc.addObservations(currentPastVectors, currentNextVectors);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.PredictiveInfoCalculator#addObservations(double[], int, int)
	 */
	@Override
	public void addObservations(double[] observations, int startTime,
			int numTimeSteps) throws Exception {
		addObservations(MatrixUtils.select(observations, startTime, numTimeSteps));
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.PredictiveInfoCalculator#finaliseAddObservations()
	 */
	@Override
	public void finaliseAddObservations() throws Exception {
		miCalc.finaliseAddObservations();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.PredictiveInfoCalculator#setObservations(double[], boolean[])
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
	 *  the embedded k-vectors with embedding delays) are valid.
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
					if (t - startTime < 2*(k-1)*tau+1) {
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
	 * @see infodynamics.measures.continuous.PredictiveInfoCalculator#computeAverageLocalOfObservations()
	 */
	@Override
	public double computeAverageLocalOfObservations() throws Exception {
		return miCalc.computeAverageLocalOfObservations();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.PredictiveInfoCalculator#computeLocalOfPreviousObservations()
	 */
	@Override
	public double[] computeLocalOfPreviousObservations() throws Exception {
		double[] local = miCalc.computeLocalOfPreviousObservations();
		if (!miCalc.getAddedMoreThanOneObservationSet()) {
			double[] localsToReturn = new double[local.length + 2*(k-1)*tau + 1];
			System.arraycopy(local, 0, localsToReturn, (k-1)*tau + 1, local.length);
			return localsToReturn;
		} else {
			return local;
		}
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.PredictiveInfoCalculator#computeLocalUsingPreviousObservations(double[])
	 */
	@Override
	public double[] computeLocalUsingPreviousObservations(double[] newObservations) throws Exception {
		if (newObservations.length - 2*(k-1)*tau - 1 <= 0) {
			// There are no observations to compute for here
			return new double[newObservations.length];
		}
		double[][] newPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(newObservations, k, tau, (k-1)*tau, newObservations.length - 2*(k-1)*tau - 1);
		double[][] newNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(newObservations, 1, 2*(k-1)*tau + 1, newObservations.length - 2*(k-1)*tau - 1);
		double[] local = miCalc.computeLocalUsingPreviousObservations(newPastVectors, newNextVectors);
		// Pad the front and back of the array with zeros where local PI isn't defined:
		double[] localsToReturn = new double[local.length + 2*(k-1)*tau + 1];
		System.arraycopy(local, 0, localsToReturn, (k-1)*tau + 1, local.length);
		return localsToReturn;

	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.PredictiveInfoCalculator#computeSignificance(int)
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		return miCalc.computeSignificance(numPermutationsToCheck);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.PredictiveInfoCalculator#computeSignificance(int[][])
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception {
		return miCalc.computeSignificance(newOrderings);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.PredictiveInfoCalculator#setDebug(boolean)
	 */
	@Override
	public void setDebug(boolean debug) {
		this.debug = debug;
		miCalc.setDebug(debug);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.PredictiveInfoCalculator#getLastAverage()
	 */
	@Override
	public double getLastAverage() {
		return miCalc.getLastAverage();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.PredictiveInfoCalculator#getNumObservations()
	 */
	@Override
	public int getNumObservations() throws Exception {
		return miCalc.getNumObservations();
	}
}
