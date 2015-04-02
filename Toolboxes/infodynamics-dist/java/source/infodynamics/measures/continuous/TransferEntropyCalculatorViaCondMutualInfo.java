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

import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;

import java.util.Vector;

/**
 * A Transfer Entropy (TE) calculator (implementing {@link TransferEntropyCalculator})
 * which is affected using a 
 * given Conditional Mutual Information (MI) calculator (implementing
 * {@link ConditionalMutualInfoCalculatorMultiVariate}) to make the calculations.
 * 
 * <p>Usage is as per the paradigm outlined for {@link TransferEntropyCalculator},
 * except that in the constructor(s) for this class the implementation for
 * a {@link ConditionalMutualInfoCalculatorMultiVariate} must be supplied.
 * </p>
 * 
 * <p>This class <i>may</i> be used directly, however users are advised that
 * several child classes are available which already plug-in the various
 * conditional MI estimators
 * to provide TE calculators (taking specific caution associated with
 * each type of estimator):</p>
 * <ul>
 * 	<li>{@link infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorGaussian}</li>
 * 	<li>{@link infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov}</li>
 * </ul>
 * 
 * TODO Delete TransferEntropyCalculatorCommon once we've switched everything over to use this?
 * Might be useful to leave it after all, and move common functionality from here to there.
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
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 */
public class TransferEntropyCalculatorViaCondMutualInfo implements
		TransferEntropyCalculator {

	/**
	 * Underlying conditional mutual information calculator
	 */
	protected ConditionalMutualInfoCalculatorMultiVariate condMiCalc;
	/**
	 * Length of past destination history to consider (embedding length)
	 */
	protected int k = 1;
	/**
	 * Embedding delay to use between elements of the destination embeding vector.
	 * We're hard-coding a delay of 1 between the history vector and the next 
	 *  observation however.
	 */
	protected int k_tau = 1;
	/**
	 * Length of past source history to consider (embedding length)
	 */
	protected int l = 1;
	/**
	 * Embedding delay to use between elements of the source embeding vector.
	 */
	protected int l_tau = 1;
	/**
	 * Source-destination next observation delay
	 */
	protected int delay = 1;
	/**
	 * Time index of the last point in the destination embedding of the first
	 *  (destination past, source past, destination next) tuple than can be 
	 *  taken from any set of time-series observations. 
	 */
	protected int startTimeForFirstDestEmbedding;

	/**
	 * Whether we're in debugging mode
	 */
	protected boolean debug = false;

	/**
	 * Construct a transfer entropy calculator using an instance of
	 * condMiCalculatorClassName as the underlying conditional mutual information calculator.
	 * 
	 * @param condMiCalculatorClassName fully qualified name of the class which must implement
	 * 	{@link ConditionalMutualInfoCalculatorMultiVariate}
	 * @throws InstantiationException if the given class cannot be instantiated
	 * @throws IllegalAccessException if illegal access occurs while trying to create an instance
	 *   of the class
	 * @throws ClassNotFoundException if the given class is not found
	 */
	public TransferEntropyCalculatorViaCondMutualInfo(String condMiCalculatorClassName)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		@SuppressWarnings("unchecked")
		Class<ConditionalMutualInfoCalculatorMultiVariate> condMiClass = 
				(Class<ConditionalMutualInfoCalculatorMultiVariate>) Class.forName(condMiCalculatorClassName);
		ConditionalMutualInfoCalculatorMultiVariate condMiCalc = condMiClass.newInstance();
		construct(condMiCalc);
	}

	/**
	 * Construct a transfer entropy calculator using an instance of
	 * condMiCalcClass as the underlying conditional mutual information calculator.
	 * 
	 * @param condMiCalcClass the class which must implement
	 * 	{@link ConditionalMutualInfoCalculatorMultiVariate}
	 * @throws InstantiationException if the given class cannot be instantiated
	 * @throws IllegalAccessException if illegal access occurs while trying to create an instance
	 *   of the class
	 * @throws ClassNotFoundException if the given class is not found
	 */
	public TransferEntropyCalculatorViaCondMutualInfo(Class<ConditionalMutualInfoCalculatorMultiVariate> condMiCalcClass)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		ConditionalMutualInfoCalculatorMultiVariate condMiCalc = condMiCalcClass.newInstance();
		construct(condMiCalc);
	}

	/**
	 * Construct this calculator by passing in a constructed but not initialised
	 * underlying Conditional Mutual information calculator.
	 * 
	 * @param condMiCalc An instantiated conditional mutual information calculator.
	 * @throws Exception if the supplied calculator has not yet been instantiated.
	 */
	public TransferEntropyCalculatorViaCondMutualInfo(ConditionalMutualInfoCalculatorMultiVariate condMiCalc) throws Exception {
		if (condMiCalc == null) {
			throw new Exception("Conditional MI calculator used to construct ConditionalTransferEntropyCalculatorViaCondMutualInfo " +
					" must have already been instantiated.");
		}
		construct(condMiCalc);
	}
	
	/**
	 * Internal method to set the conditional mutual information calculator.
	 * Can be overridden if anything else needs to be done with it by the child classes.
	 * 
	 * @param condMiCalc
	 */
	protected void construct(ConditionalMutualInfoCalculatorMultiVariate condMiCalc) {
		this.condMiCalc = condMiCalc;
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ChannelCalculatorCommon#initialise()
	 */
	@Override
	public void initialise() throws Exception {
		initialise(k, k_tau, l, l_tau, delay);
	}
	
	@Override
	public void initialise(int k) throws Exception {
		initialise(k, k_tau, l, l_tau, delay);
	}
	
	/**
	 * Initialise the calculator for re-use with new observations.
	 * New embedding parameters and source-destination delay
	 * may be supplied here; all other parameters
	 * remain unchanged.
	 * 
	 * @param k embedding length of destination past history to consider
	 * @param k_tau embedding delay for the destination variable
	 * @param l embedding length of source past history to consider
	 * @param l_tau embedding delay for the source variable
	 * @param delay time lag between last element of source and destination next value
	 */
	public void initialise(int k, int k_tau, int l, int l_tau, int delay) throws Exception {
		if (delay < 0) {
			throw new Exception("Cannot compute TE with source-destination delay < 0");
		}
		this.k = k;
		this.k_tau = k_tau;
		this.l = l;
		this.l_tau = l_tau;
		this.delay = delay;
		
		// Now check which point we can start taking observations from in any
		//  addObservations call. These two integers represent the last
		//  point of the destination embedding, in the cases where the destination
		//  embedding itself determines where we can start taking observations, or
		//  the case where the source embedding plus delay is longer and so determines
		//  where we can start taking observations.
		int startTimeBasedOnDestPast = (k-1)*k_tau;
		int startTimeBasedOnSourcePast = (l-1)*l_tau + delay - 1;
		startTimeForFirstDestEmbedding = Math.max(startTimeBasedOnDestPast, startTimeBasedOnSourcePast);

		condMiCalc.initialise(l, 1, k);
	}

	/**
	 * Sets properties for the TE calculator.
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 *  
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 * 		<li>Any properties accepted by {@link TransferEntropyCalculator#setProperty(String, String)}</li>
	 * 		<li>Or properties accepted by the underlying
	 * 		{@link ConditionalMutualInfoCalculatorMultiVariate#setProperty(String, String)} implementation.</li>
	 * </ul>
	 * <p><b>Note:</b> further properties may be defined by child classes.</p>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property.
	 * @throws Exception if there is a problem with the supplied value.
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(K_PROP_NAME)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(K_TAU_PROP_NAME)) {
			k_tau = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(L_PROP_NAME)) {
			l = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(L_TAU_PROP_NAME)) {
			l_tau = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(DELAY_PROP_NAME)) {
			delay = Integer.parseInt(propertyValue);
		} else {
			// No property was set on this class, assume it is for the underlying
			//  conditional MI calculator
			condMiCalc.setProperty(propertyName, propertyValue);
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	@Override
	public void setObservations(double[] source, double[] destination) throws Exception {
		
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		if (source.length < startTimeForFirstDestEmbedding + 2) {
			// There are no observations to add here, the time series is too short
			throw new Exception("Not enough observations to set here given k, k_tau, l, l_tau and delay parameters");
		}
		double[][] currentDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(destination, k, k_tau,
						startTimeForFirstDestEmbedding,
						destination.length - startTimeForFirstDestEmbedding - 1);
		double[][] currentDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(destination, 1,
						startTimeForFirstDestEmbedding + 1,
						destination.length - startTimeForFirstDestEmbedding - 1);
		double[][] currentSourcePastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(source, l, l_tau,
						startTimeForFirstDestEmbedding + 1 - delay,
						source.length - startTimeForFirstDestEmbedding - 1);
		condMiCalc.setObservations(currentSourcePastVectors, currentDestNextVectors, currentDestPastVectors);
	}

	@Override
	public void startAddObservations() {
		condMiCalc.startAddObservations();
	}
	
	@Override
	public void addObservations(double[] source, double[] destination) throws Exception {
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		if (source.length < startTimeForFirstDestEmbedding + 2) {
			// There are no observations to add here, the time series is too short
			// Don't throw an exception, do nothing since more observations
			//  can be added later.
			return;
		}
		double[][] currentDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(destination, k, k_tau,
						startTimeForFirstDestEmbedding,
						destination.length - startTimeForFirstDestEmbedding - 1);
		double[][] currentDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(destination, 1,
						startTimeForFirstDestEmbedding + 1,
						destination.length - startTimeForFirstDestEmbedding - 1);
		double[][] currentSourcePastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(source, l, l_tau,
						startTimeForFirstDestEmbedding + 1 - delay,
						source.length - startTimeForFirstDestEmbedding - 1);
		condMiCalc.addObservations(currentSourcePastVectors, currentDestNextVectors, currentDestPastVectors);
	}

	@Override
	public void addObservations(double[] source, double[] destination,
			int startTime, int numTimeSteps) throws Exception {
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		if (source.length < startTime + numTimeSteps) {
			// There are not enough observations given the arguments here
			throw new Exception("Not enough observations to set here given startTime and numTimeSteps parameters");
		}
		addObservations(MatrixUtils.select(source, startTime, numTimeSteps),
					    MatrixUtils.select(destination, startTime, numTimeSteps));
	}

	@Override
	public void finaliseAddObservations() throws Exception {
		condMiCalc.finaliseAddObservations();
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
	 *  valid series of both source and destinations. (I.e. all points within the
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
	public Vector<int[]> computeStartAndEndTimePairs(boolean[] sourceValid, boolean[] destValid) throws Exception {
		
		if (sourceValid.length != destValid.length) {
			throw new Exception("Validity arrays must be of same length");
		}
		
		int lengthOfDestPastRequired = (k-1)*k_tau + 1;
		int lengthOfSourcePastRequired = (l-1)*l_tau + 1;
		// int numSourcePointsBeforeDestStart = delay - 1 + lengthOfSourcePastRequired
		//									- lengthOfDestPastRequired;

		// Scan along the data avoiding invalid values
		int startTime = 0;
		Vector<int[]> startAndEndTimePairs = new Vector<int[]>();

		// Simple solution -- this takes more complexity in time, but is 
		//  much faster to code:
		boolean previousWasOk = false;
		for (int t = startTimeForFirstDestEmbedding; t < destValid.length - 1; t++) {
			// Check the tuple with the history vector starting from
			//  t and running backwards
			if (previousWasOk) {
				// Just check the very next values of each:
				if (destValid[t + 1] && sourceValid[t + 1 - delay]) {
					// We can continue adding to this sequence
					continue;
				} else {
					// We need to shut down this sequence now
					previousWasOk = false;
					int[] timePair = new int[2];
					timePair[0] = startTime;
					timePair[1] = t; // Previous time step was last valid one
					startAndEndTimePairs.add(timePair);
					continue;
				}
			}
			// Otherwise we're trying to start a new sequence, so check all values
			if (!destValid[t + 1]) {
				continue;
			}
			boolean allOk = true;
			for (int tBack = 0; tBack < lengthOfDestPastRequired; tBack++) {
				if (!destValid[t - tBack]) {
					allOk = false;
					break;
				}
			}
			if (!allOk) {
				continue;
			}
			allOk = true;
			for (int tBack = delay - 1; tBack < delay - 1 + lengthOfSourcePastRequired; tBack++) {
				if (!sourceValid[t - tBack]) {
					allOk = false;
					break;
				}
			}
			if (!allOk) {
				continue;
			}
			// Postcondition: We've got a first valid tuple:
			startTime = t - startTimeForFirstDestEmbedding;
			previousWasOk = true;
		}
		// Now check if we were running a sequence and terminate it:
		if (previousWasOk) {
			// We need to shut down this sequence now
			previousWasOk = false;
			int[] timePair = new int[2];
			timePair[0] = startTime;
			timePair[1] = destValid.length - 1;
			startAndEndTimePairs.add(timePair);
		}

		return startAndEndTimePairs;
	}
	
	@Override
	public double computeAverageLocalOfObservations() throws Exception {
		return condMiCalc.computeAverageLocalOfObservations();
	}

	@Override
	public double[] computeLocalOfPreviousObservations() throws Exception {
		double[] local = condMiCalc.computeLocalOfPreviousObservations();
		if (!condMiCalc.getAddedMoreThanOneObservationSet()) {
			double[] localsToReturn = new double[local.length + startTimeForFirstDestEmbedding + 1];
			System.arraycopy(local, 0, localsToReturn, startTimeForFirstDestEmbedding + 1, local.length);
			return localsToReturn;
		} else {
			return local;
		}
	}

	@Override
	public double[] computeLocalUsingPreviousObservations(
			double[] newSourceObservations, double[] newDestObservations)
					throws Exception {
		if (newSourceObservations.length != newDestObservations.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					newSourceObservations.length, newDestObservations.length));
		}
		if (newDestObservations.length < startTimeForFirstDestEmbedding + 2) {
			// There are no observations to compute for here
			return new double[newDestObservations.length];
		}
		double[][] newDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(newDestObservations, k, k_tau,
						startTimeForFirstDestEmbedding,
						newDestObservations.length - startTimeForFirstDestEmbedding - 1);
		double[][] newDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(newDestObservations, 1,
						startTimeForFirstDestEmbedding + 1,
						newDestObservations.length - startTimeForFirstDestEmbedding - 1);
		double[][] newSourcePastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(newSourceObservations, l, l_tau,
						startTimeForFirstDestEmbedding + 1 - delay,
						newSourceObservations.length - startTimeForFirstDestEmbedding - 1);
		double[] local = condMiCalc.computeLocalUsingPreviousObservations(
						newSourcePastVectors, newDestNextVectors, newDestPastVectors);
		// Pad the front of the array with zeros where local TE isn't defined:
		double[] localsToReturn = new double[local.length + startTimeForFirstDestEmbedding + 1];
		System.arraycopy(local, 0, localsToReturn, startTimeForFirstDestEmbedding + 1, local.length);
		return localsToReturn;

	}
	
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		// Reorder the source vectors in the surrogates, not the destination
		return condMiCalc.computeSignificance(1, numPermutationsToCheck);
	}

	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception {
		// Reorder the source vectors in the surrogates, not the destination
		return condMiCalc.computeSignificance(1, newOrderings);
	}

	@Override
	public double getLastAverage() {
		return condMiCalc.getLastAverage();
	}

	@Override
	public int getNumObservations() throws Exception {
		return condMiCalc.getNumObservations();
	}
	
	@Override
	public boolean getAddedMoreThanOneObservationSet() {
		return condMiCalc.getAddedMoreThanOneObservationSet();
	}

	@Override
	public void setDebug(boolean debug) {
		this.debug = debug;
		condMiCalc.setDebug(debug);
	}
}
