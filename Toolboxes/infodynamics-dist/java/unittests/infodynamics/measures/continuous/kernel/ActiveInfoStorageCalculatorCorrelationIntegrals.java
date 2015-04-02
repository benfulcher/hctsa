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
import infodynamics.utils.MatrixUtils;

import java.util.Vector;


public abstract class ActiveInfoStorageCalculatorCorrelationIntegrals
	implements ActiveInfoStorageCalculator {

	/**
	 * Length of past history to consider
	 */
	protected int k = 1;
	protected int totalObservations = 0;
	protected boolean debug = false;
	protected double lastAverage;
	
	/**
	 * Storage for source observations for addObservsations
	 */
	protected Vector<double[]> vectorOfObservations;
	
	protected boolean addedMoreThanOneObservationSet;

	public ActiveInfoStorageCalculatorCorrelationIntegrals() {
	}

	/**
	 * Initialise the calculator using the existing value of k
	 * 
	 */
	public void initialise() throws Exception {
		initialise(k);
	}

	/**
	 * Initialise the calculator
	 * 
	 * @param k Length of past history to consider
	 */
	public void initialise(int k) throws Exception {
		this.k = k;
		addedMoreThanOneObservationSet = false;
	}
	
	/**
	 * Set the observations to compute the probabilities from 
	 * 
	 * @param observations
	 */
	public void setObservations(double[] observations) throws Exception {
		startAddObservations();
		addObservations(observations);
		finaliseAddObservations();
	}

	/**
	 * Elect to add in the observations from several disjoint time series.
	 *
	 */
	public void startAddObservations() {
		vectorOfObservations = new Vector<double[]>();
	}
	
	/**
	 * Add some more observations.
	 * Note that the arrays must not be over-written by the user
	 *  until after finaliseAddObservations() has been called.
	 * 
	 * @param observations
	 */
	public void addObservations(double[] observations) throws Exception {
		if (vectorOfObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		if (observations.length <= k) {
			// we won't be taking any observations here
			return;
		}
		vectorOfObservations.add(observations);
	}

	/**
	 * Add some more observations.
	 * 
	 * @param observations
	 * @param startTime first time index to take observations on
	 * @param numTimeSteps number of time steps to use
	 */
	public void addObservations(double[] observations,
			int startTime, int numTimeSteps) throws Exception {
		if (vectorOfObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		if (numTimeSteps <= k) {
			// We won't be taking any observations here
			return;
		}
		double[] obsToAdd = new double[numTimeSteps];
		System.arraycopy(observations, startTime, obsToAdd, 0, numTimeSteps);
		vectorOfObservations.add(obsToAdd);
	}

	/**
	 * Sets the observations to compute the PDFs from.
	 * Cannot be called in conjunction with start/add/finaliseAddObservations.
	 * valid is a time series (with time indices the same as observations)
	 *  indicating whether the observation at that point is valid.
	 * sourceValid is the same for the source
	 * 
	 * @param observation observations for the source variable
	 * @param valid
	 */
	public void setObservations(double[] observations,
			boolean[] valid) throws Exception {
		
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
	 *  valid series of the observations.
	 * 
	 * Made public so it can be used if one wants to compute the number of
	 *  observations prior to setting the observations.
	 * 
	 * @param valid
	 * @return
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
					if (t - startTime < k) {
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
	
	/**
	 * Generate a vector for each time step, containing the past k states of the destination.
	 * Does not include a vector for the first k time steps.
	 * 
	 * @param destination
	 * @return array of vectors for each time step
	 */
	protected double[][] makeJointVectorForPast(double[] destination) {
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
	
	/**
	 * Compute the next values into a vector
	 * 
	 * @param destination
	 * @return
	 */
	protected double[][] makeJointVectorForNext(double[] destination) {
		double[][] destNextVectors = new double[destination.length - k][1];
		for (int t = k; t < destination.length; t++) {
			destNextVectors[t - k][0] = destination[t];
		}
		return destNextVectors;
	}
	
	/**
	 * Generate a vector for each time step, containing the past k states of
	 *  the destination, and the current state.
	 * Does not include a vector for the first k time steps.
	 * 
	 * @param destination
	 * @return
	 */
	protected double[][] makeJointVectorForNextPast(double[] destination) {
		// We want all delay vectors here
		return MatrixUtils.makeDelayEmbeddingVector(destination, k+1);
	}
	

	public double getLastAverage() {
		return lastAverage;
	}

	public int getNumObservations() {
		return totalObservations;
	}
	
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(K_PROP_NAME)) {
			k = Integer.parseInt(propertyValue);
		} else {
			// No property was set
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}
}
