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

import infodynamics.utils.MatrixUtils;

import java.util.Vector;

/**
 * Implements {@link TransferEntropyCalculator} to provide a base
 * class with common functionality for child class implementations of
 * pairwise or apparent transfer entropy (TE), defined by {@link TransferEntropyCalculator},
 * on <code>double[]</code> data via various estimators. 
 * 
 * <p>These various estimators include: e.g. box-kernel estimation, KSG estimators, etc
 * (see the child classes linked above).
 * </p>
 * 
 * <p>Usage is as outlined in {@link TransferEntropyCalculator}.</p>
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
public abstract class TransferEntropyCommon implements
		TransferEntropyCalculator {

	/**
	 * Natural log of 2, cached
	 */
	protected final static double log2 = Math.log(2.0);
	
	/**
	 * Length of past destination history to consider (embedding length)
	 */
	protected int k;
	/**
	 * Total number of observations supplied.
	 * Only valid after {@link #finaliseAddObservations()} is called.
	 */
	protected int totalObservations = 0;
	/**
	 * Whether to report debug messages or not
	 */
	protected boolean debug = false;
	/**
	 * Store the last computed average TE
	 */
	protected double lastAverage;

	/**
	 * Storage for source observations for {@link #addObservations(double[], double[])} etc
	 */
	protected Vector<double[]> vectorOfSourceObservations;
	/**
	 * Storage for destination observations for {@link #addObservations(double[], double[])} etc
	 */
	protected Vector<double[]> vectorOfDestinationObservations;
	
	/**
	 * Whether the user has supplied more than one (disjoint) set of samples
	 */
	protected boolean addedMoreThanOneObservationSet;

	@Override
	public void initialise(int k) throws Exception {
		this.k = k;
		addedMoreThanOneObservationSet = false;
		initialise();
	}

	@Override
	public void setObservations(double[] source, double[] destination) throws Exception {
		startAddObservations();
		addObservations(source, destination);
		finaliseAddObservations();
	}

	@Override
	public void startAddObservations() {
		vectorOfSourceObservations = new Vector<double[]>();
		vectorOfDestinationObservations = new Vector<double[]>();
	}
	
	@Override
	public void addObservations(double[] source, double[] destination) throws Exception {
		if (vectorOfSourceObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		if (source.length <= k) {
			// we won't be taking any observations here
			return;
		}
		vectorOfSourceObservations.add(source);
		vectorOfDestinationObservations.add(destination);
	}

	@Override
	public void addObservations(double[] source, double[] destination,
			int startTime, int numTimeSteps) throws Exception {
		if (vectorOfSourceObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		if (numTimeSteps <= k) {
			// We won't be taking any observations here
			return;
		}
		double[] sourceToAdd = new double[numTimeSteps];
		System.arraycopy(source, startTime, sourceToAdd, 0, numTimeSteps);
		vectorOfSourceObservations.add(sourceToAdd);
		double[] destToAdd = new double[numTimeSteps];
		System.arraycopy(destination, startTime, destToAdd, 0, numTimeSteps);
		vectorOfDestinationObservations.add(destToAdd);
	}

	@Override
	public void setObservations(double[] source, double[] destination,
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

	/**
	 * Compute a vector of start and end pairs of time points, between which we have
	 *  valid series of both source and destinations.  (I.e. all points within the
	 *  embedding vectors must be valid, even if the invalid points won't be included
	 *  in any tuples)
	 * 
	 * <p>Made public so it can be used if one wants to compute the number of
	 *  observations prior to setting the observations.</p>
	 * 
	 * @param sourceValid a time series (with indices the same as observations)
	 *  indicating whether the entry in observations at that index is valid for the source; 
	 * @param destValid as described for <code>sourceValid</code>
	 * @return a vector for start and end time pairs of valid series
	 *  of observations.
	 */
	public Vector<int[]> computeStartAndEndTimePairs(boolean[] sourceValid, boolean[] destValid) {
		// Scan along the data avoiding invalid values
		int startTime = 0;
		int endTime = 0;
		boolean lookingForStart = true;
		Vector<int[]> startAndEndTimePairs = new Vector<int[]>();
		for (int t = 0; t < destValid.length; t++) {
			if (lookingForStart) {
				// Precondition: startTime holds a candidate start time
				if (destValid[t]) {
					// This point is OK at the destination
					if (t - startTime < k) {
						// We're still checking the past history only, so
						continue;
					} else {
						// We've got the full past history ok, so check the source also
						if (sourceValid[t - 1]) {
							// source is good to go also
							// set a candidate endTime
							endTime = t;
							lookingForStart = false;
							if (t == destValid.length - 1) {
								// we need to terminate now
								int[] timePair = new int[2];
								timePair[0] = startTime;
								timePair[1] = endTime;
								startAndEndTimePairs.add(timePair);
								// System.out.printf("t_s=%d, t_e=%d\n", startTime, endTime);
							}
						} else {
							// source was not good to go, so try moving along one time point
							startTime++;
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
				if (destValid[t] && sourceValid[t - 1]) {
					// We can extend
					endTime = t;
				} else {
					terminateSequence = true;
				}
				if (t == destValid.length - 1) {
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
					if (!destValid[t]) {
						// The current destination observation broke our chain;
						//  so we need to start looking all over again:
						startTime = t + 1;
					} else {
						// The current source broke our chain (or we're at the
						//  end of the series anyway, so this doesn't matter);
						//  so we can keep the good destination history
						//  that we've built up here:
						startTime = t - k + 1;
					}
				}
			}
		}
		return startAndEndTimePairs;
	}
	
	/**
	 * Generate an embedding vector for each time step, containing the past k states of the destination.
	 * Does not include a vector for the first k time steps.
	 * 
	 * @param destination time series for the destination variable
	 * @return array of embedding vectors for each time step (of length destination.length - k)
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
	 * Generate an embedding vector for each time step, containing the past k+1 states of
	 *  the destination
	 * Does not include a vector for the first k time steps.
	 * 
	 * @param destination time series for the destination variable
	 * @return array of embedding vectors for each time step (of length destination.length - k)
	 */
	protected double[][] makeJointVectorForNextPast(double[] destination) {
		// We want all delay vectors here
		return MatrixUtils.makeDelayEmbeddingVector(destination, k+1);
	}
	
	/**
	 * Sets properties for the TE calculator.
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 *  
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 * 		<li>{@link #K_PROP_NAME} -- specified by {@link TransferEntropyCalculator}</li>
	 * </ul>
	 * <p><b>However -- </b> any other properties specified by {@link TransferEntropyCalculator}
	 *  (i.e. {@link #K_TAU_PROP_NAME}, {@link #L_PROP_NAME}, {@link #L_TAU_PROP_NAME} 
	 *  or {@link #DELAY_PROP_NAME}) are not implemented and will cause 
	 *  an exception to be thrown.</p>
	 *  
	 * <p>Unknown property values are otherwise ignored.</p>
	 * 
	 * <p><b>Note:</b> further properties may be defined by child classes.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property.
	 * @throws Exception if there is a problem with the supplied value, 
	 * or if the property is recognised but unsupported (as described above).
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(K_PROP_NAME)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(K_TAU_PROP_NAME)) {
			throw new Exception("Unsupported property");
		} else if (propertyName.equalsIgnoreCase(L_PROP_NAME)) {
			throw new Exception("Unsupported property");
		} else if (propertyName.equalsIgnoreCase(L_TAU_PROP_NAME)) {
			throw new Exception("Unsupported property");
		} else if (propertyName.equalsIgnoreCase(DELAY_PROP_NAME)) {
			throw new Exception("Unsupported property");
		} else {
			// No property was set
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	@Override
	public double getLastAverage() {
		return lastAverage;
	}

	@Override
	public int getNumObservations() {
		return totalObservations;
	}
	
	@Override
	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	@Override
	public boolean getAddedMoreThanOneObservationSet() {
		return addedMoreThanOneObservationSet;
	}
}
