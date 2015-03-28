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
import infodynamics.utils.RandomGenerator;

import java.util.Iterator;
import java.util.Vector;

/**
 * Implements {@link ConditionalMutualInfoCalculatorMultiVariate}
 * to provide a base
 * class with common functionality for child class implementations of
 * {@link ConditionalMutualInfoCalculatorMultiVariate}
 * via various estimators. 
 * 
 * <p>These various estimators include: e.g. box-kernel estimation, KSG estimators, etc
 * (see the child classes linked above).
 * </p>
 * 
 * <p>Usage is as outlined in {@link ConditionalMutualInfoCalculatorMultiVariate}.</p>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public abstract class ConditionalMutualInfoMultiVariateCommon implements
		ConditionalMutualInfoCalculatorMultiVariate {

	/**
	 * Number of dimenions for variable 1
	 */
	protected int dimensionsVar1 = 1;
	/**
	 * Number of dimenions for variable 2
	 */
	protected int dimensionsVar2 = 1;
	/**
	 * Number of dimenions for the conditional variable
	 */
	protected int dimensionsCond = 1;
	
	/**
	 * The set of observations for var1, retained in case the user wants to retrieve the local
	 *  entropy values of these.
	 * They're held in the order in which they were supplied in the
	 *  {@link #addObservations(double[][], double[][], double[][])} functions.
	 */
	protected double[][] var1Observations;

	/**
	 * The set of observations for var2, retained in case the user wants to retrieve the local
	 *  entropy values of these.
	 * They're held in the order in which they were supplied in the
	 *  {@link #addObservations(double[][], double[][], double[][])} functions.
	 */
	protected double[][] var2Observations;

	/**
	 * The set of observations for the conditional, retained in case the user wants to retrieve the local
	 *  entropy values of these.
	 * They're held in the order in which they were supplied in the
	 *  {@link #addObservations(double[][], double[][], double[][])} functions.
	 */
	protected double[][] condObservations;

	/**
	 * Total number of observations supplied.
	 * Only valid after {@link #finaliseAddObservations()} is called.
	 */
	protected int totalObservations = 0;

	/**
	 * Store the last computed average conditional MI
	 */
	protected double lastAverage;

	/**
	 * Track whether we've computed the average for the supplied
	 *  observations yet
	 */
	protected boolean condMiComputed;

	/**
	 * Whether to report debug messages or not
	 */
	protected boolean debug;

	/**
	 * Storage for var1 observations supplied via
	 * {@link #addObservations(double[][], double[][], double[][])} etc
	 */
	protected Vector<double[][]> vectorOfVar1Observations;
	/**
	 * Storage for var2 observations supplied via
	 * {@link #addObservations(double[][], double[][], double[][])} etc
	 */
	protected Vector<double[][]> vectorOfVar2Observations;
	/**
	 * Storage for conditional variable observations supplied via
	 * {@link #addObservations(double[][], double[][], double[][])} etc
	 */
	protected Vector<double[][]> vectorOfCondObservations;
	
	/**
	 * Whether the user has added more than one disjoint observation set
	 * via {@link #addObservations(double[][], double[][], double[][])} etc
	 */
	protected boolean addedMoreThanOneObservationSet;

	/**
	 * Property name for 
	 *  specifying whether the data is normalised or not (to mean 0,
	 *  variance 1, for each of the multiple variables)
	 *  before the calculation is made.
	 */
	public static final String PROP_NORMALISE = "NORMALISE";
	protected boolean normalise = true;

	@Override
	public void initialise(int var1Dimensions, int var2Dimensions, int condDimensions) {
		dimensionsVar1 = var1Dimensions;
		dimensionsVar2 = var2Dimensions;
		dimensionsCond = condDimensions;
		lastAverage = 0.0;
		totalObservations = 0;
		condMiComputed = false;
		var1Observations = null;
		var2Observations = null;
		condObservations = null;
		addedMoreThanOneObservationSet = false;
	}

	/**
	 * Set properties for the calculator.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 *  	<li>{@link #PROP_NORMALISE} - whether to normalise the individual
	 *      variables to mean 0, standard deviation 1
	 *      (true by default, except for child class
	 *      {@link infodynamics.measures.continuous.gaussian.ConditionalMutualInfoCalculatorMultiVariateGaussian})</li>
	 * </ul>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue) {
		if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
			normalise = Boolean.parseBoolean(propertyValue);
		}
	}

	@Override
	public void setObservations(double[][] var1, double[][] var2,
			double[][] cond) throws Exception {
		startAddObservations();
		addObservations(var1, var2, cond);
		finaliseAddObservations();
		addedMoreThanOneObservationSet = false;
	}

	@Override
	public void startAddObservations() {
		vectorOfVar1Observations = new Vector<double[][]>();
		vectorOfVar2Observations = new Vector<double[][]>();
		vectorOfCondObservations = new Vector<double[][]>();
	}
	
	@Override
	public void addObservations(double[][] var1, double[][] var2,
			double[][] cond) throws Exception {
		if (vectorOfVar1Observations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		if ((var1.length != var2.length) || (var1.length != cond.length)) {
			throw new Exception(String.format("Observation vector lengths (%d, %d and %d) must match!",
					var1.length, var2.length, cond.length));
		}
		if (var1[0].length != dimensionsVar1) {
			throw new Exception("Number of joint variables in var1 data " +
					"does not match the initialised value");
		}
		if (var2[0].length != dimensionsVar2) {
			throw new Exception("Number of joint variables in var2 data " +
					"does not match the initialised value");
		}
		if (cond[0].length != dimensionsCond) {
			throw new Exception("Number of joint variables in cond data " +
					"does not match the initialised value");
		}
		vectorOfVar1Observations.add(var1);
		vectorOfVar2Observations.add(var2);
		vectorOfCondObservations.add(cond);
		if (vectorOfVar1Observations.size() > 1) {
			addedMoreThanOneObservationSet = true;
		}
	}

	@Override
	public void addObservations(double[][] var1, double[][] var2,
			double[][] cond,
			int startTime, int numTimeSteps) throws Exception {
		if (vectorOfVar1Observations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		double[][] var1ToAdd = new double[numTimeSteps][];
		System.arraycopy(var1, startTime, var1ToAdd, 0, numTimeSteps);
		double[][] var2ToAdd = new double[numTimeSteps][];
		System.arraycopy(var2, startTime, var2ToAdd, 0, numTimeSteps);
		double[][] condToAdd = new double[numTimeSteps][];
		System.arraycopy(cond, startTime, condToAdd, 0, numTimeSteps);
		addObservations(var1ToAdd, var2ToAdd, condToAdd);
	}

	@Override
	public void setObservations(double[][] var1, double[][] var2,
			double[][] cond,
			boolean[] var1Valid, boolean[] var2Valid,
			boolean[] condValid) throws Exception {
		
		Vector<int[]> startAndEndTimePairs =
				computeStartAndEndTimePairs(var1Valid, var2Valid, condValid);
		
		// We've found the set of start and end times for this pair
		startAddObservations();
		for (int[] timePair : startAndEndTimePairs) {
			int startTime = timePair[0];
			int endTime = timePair[1];
			addObservations(var1, var2, cond, startTime, endTime - startTime + 1);
		}
		finaliseAddObservations();
	}

	@Override
	public void setObservations(double[][] var1, double[][] var2,
			double[][] cond,
			boolean[][] var1Valid, boolean[][] var2Valid,
			boolean[][] condValid) throws Exception {

		boolean[] allVar1Valid = MatrixUtils.andRows(var1Valid);
		boolean[] allVar2Valid = MatrixUtils.andRows(var2Valid);
		boolean[] allCondValid = MatrixUtils.andRows(condValid);
		setObservations(var1, var2, cond, allVar1Valid, allVar2Valid, allCondValid);
	}

	/**
	 * Signal that the observations are now all added, PDFs can now be constructed.
	 * 
	 * This default implementation simply puts all of the observations into
	 *  the {@link #var1Observations}, {@link #var2Observations} 
	 *  and {@link #condObservations} arrays.
	 * Usually child implementations will override this, call this implementation
	 *  to perform the common processing, then perform their own processing.
	 *  
	 * @throws Exception Allow child classes to throw an exception if there
	 *  is an issue detected specific to that calculator.
	 */
	@Override
	public void finaliseAddObservations() throws Exception {
		// First work out the size to allocate the joint vectors, and do the allocation:
		totalObservations = 0;
		for (double[][] var2 : vectorOfVar2Observations) {
			totalObservations += var2.length;
		}
		var1Observations = new double[totalObservations][dimensionsVar1];
		var2Observations = new double[totalObservations][dimensionsVar2];
		condObservations = new double[totalObservations][dimensionsCond];
		
		int startObservation = 0;
		Iterator<double[][]> iteratorVar2 = vectorOfVar2Observations.iterator();
		Iterator<double[][]> iteratorCond = vectorOfCondObservations.iterator();
		for (double[][] var1 : vectorOfVar1Observations) {
			double[][] var2 = iteratorVar2.next();
			double[][] cond = iteratorCond.next();
			// Copy the data from these given observations into our master 
			//  array
			MatrixUtils.arrayCopy(var1, 0, 0,
					var1Observations, startObservation, 0,
					var1.length, dimensionsVar1);
			MatrixUtils.arrayCopy(var2, 0, 0,
					var2Observations, startObservation, 0,
					var2.length, dimensionsVar2);
			MatrixUtils.arrayCopy(cond, 0, 0,
					condObservations, startObservation, 0,
					cond.length, dimensionsCond);
			startObservation += var2.length;
		}
		
		// Normalise the data if required
		if (normalise) {
			MatrixUtils.normalise(var1Observations);
			MatrixUtils.normalise(var2Observations);
			MatrixUtils.normalise(condObservations);
		}
		
		// We don't need to keep the vectors of observation sets anymore:
		vectorOfVar1Observations = null;
		vectorOfVar2Observations = null;
		vectorOfCondObservations = null;
	}
	
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int variableToReorder, int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		// Use var1 length (all variables have same length) even though
		//  we may be randomising the other variable:
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(
				var1Observations.length, numPermutationsToCheck);
		return computeSignificance(variableToReorder, newOrderings);
	}

	/**
	 * <p>As described in
	 * {@link ConditionalMutualInfoCalculatorMultiVariate#computeSignificance(int, int[][])}
	 * </p>
	 * 
	 * <p>Here we provide a simple implementation which would be suitable for
	 *  any child class, though the child class may prefer to make its
	 *  own implementation to make class-specific optimisations.
	 *  Child classes must implement {@link java.lang.Cloneable}
	 *  for this method to be callable for them, and indeed implement
	 *  the clone() method in a way that protects their structure
	 *  from alteration by surrogate data being supplied to it.</p>
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int variableToReorder, int[][] newOrderings) throws Exception {
		
		int numPermutationsToCheck = newOrderings.length;
		if (!condMiComputed) {
			computeAverageLocalOfObservations();
		}
		
		// Take a clone of the object to compute the MI of the surrogates:
		// (this is a shallow copy, it doesn't make new copies of all
		//  the arrays - child classes should override this)
		ConditionalMutualInfoMultiVariateCommon miSurrogateCalculator =
				(ConditionalMutualInfoMultiVariateCommon) this.clone();
		
		double[] surrogateMeasurements = new double[numPermutationsToCheck];
		
		// Now compute the MI for each set of shuffled data:
		for (int i = 0; i < numPermutationsToCheck; i++) {
			// Generate a new re-ordered source data
			double[][] shuffledData = 
					MatrixUtils.extractSelectedTimePointsReusingArrays(
							(variableToReorder == 1) ? var1Observations : var2Observations,
							newOrderings[i]);
			// Perform new initialisations
			miSurrogateCalculator.initialise(
					dimensionsVar1, dimensionsVar2, dimensionsCond);
			// Set new observations
			if (variableToReorder == 1) {
				miSurrogateCalculator.setObservations(shuffledData,
						var2Observations, condObservations);
			} else {
				miSurrogateCalculator.setObservations(var1Observations,
						shuffledData, condObservations);
			}
			// Compute the MI
			surrogateMeasurements[i] = miSurrogateCalculator.computeAverageLocalOfObservations();
			if (debug){
				System.out.println("New MI was " + surrogateMeasurements[i]);
			}
		}
		
		return new EmpiricalMeasurementDistribution(surrogateMeasurements, lastAverage);
	}

	/**
	 * <p>As described in
	 * {@link ConditionalMutualInfoCalculatorMultiVariate#computeAverageLocalOfObservations(int, int[])}
	 * </p>
	 * 
	 * <p>We provide a simple implementation which would be suitable for
	 *  any child class, though the child class may prefer to make its
	 *  own implementation to make class-specific optimisations.
	 *  Child classes must implement {@link java.lang.Cloneable}
	 *  for this method to be callable for them, and indeed implement
	 *  the clone() method in a way that protects their structure
	 *  from alteration by surrogate data being supplied to it.</p>
	 */
	@Override
	public double computeAverageLocalOfObservations(int variableToReorder, int[] newOrdering)
			throws Exception {
		// Take a clone of the object to compute the MI of the surrogates:
		// (this is a shallow copy, it doesn't make new copies of all
		//  the arrays - child class should override this)
		ConditionalMutualInfoMultiVariateCommon miSurrogateCalculator =
				(ConditionalMutualInfoMultiVariateCommon) this.clone();
		
		// Generate a new re-ordered source data
		double[][] shuffledData =
				MatrixUtils.extractSelectedTimePointsReusingArrays(
					(variableToReorder == 1) ? var1Observations : var2Observations,
					newOrdering);
		// Perform new initialisations
		miSurrogateCalculator.initialise(
				dimensionsVar1, dimensionsVar2, dimensionsCond);
		// Set new observations
		if (variableToReorder == 1) {
			miSurrogateCalculator.setObservations(shuffledData,
					var2Observations, condObservations);
		} else {
			miSurrogateCalculator.setObservations(var1Observations,
					shuffledData, condObservations);
		}
		// Compute the MI
		return miSurrogateCalculator.computeAverageLocalOfObservations();
	}

	@Override
	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public double getLastAverage() {
		return lastAverage;
	}

	public int getNumObservations() throws Exception {
		return totalObservations;
	}

	/**
	 * Compute a vector of start and end pairs of time points, between which we have
	 *  valid series of all variables.
	 * 
	 * Made public so it can be used if one wants to compute the number of
	 *  observations prior to setting the observations.
	 * 
	 * @param var1Valid a series (indexed by observation number or time)
	 *  indicating whether the entry in observations at that index is valid for variable 1; 
	 * @param var2Valid as described for <code>var1Valid</code>
	 * @param condValid as described for <code>var1Valid</code>
	 * @return a vector for start and end time pairs of valid series
	 *  of observations.
	 */
	public Vector<int[]> computeStartAndEndTimePairs(
			boolean[] var1Valid, boolean[] var2Valid, boolean[] condValid) {
		// Scan along the data avoiding invalid values
		int startTime = 0;
		int endTime = 0;
		boolean lookingForStart = true;
		Vector<int[]> startAndEndTimePairs = new Vector<int[]>();
		for (int t = 0; t < var2Valid.length; t++) {
			if (lookingForStart) {
				// Precondition: startTime holds a candidate start time
				//  (var1 value is at startTime == t)
				if (var1Valid[t] && var2Valid[t] && condValid[t]) {
					// This point is OK at the variables
					// Set a candidate endTime
					endTime = t;
					lookingForStart = false;
					if (t == var1Valid.length - 1) {
						// we need to terminate now
						int[] timePair = new int[2];
						timePair[0] = startTime;
						timePair[1] = endTime;
						startAndEndTimePairs.add(timePair);
						// System.out.printf("t_s=%d, t_e=%d\n", startTime, endTime);
					}
				} else {
					// We need to keep looking.
					// Move the potential start time to the next point
					startTime++;
				}
			} else {
				// Precondition: startTime holds the start time for this set, 
				//  endTime holds a candidate end time
				// Check if we can include the current time step
				boolean terminateSequence = false;
				if (var1Valid[t] && var2Valid[t] && condValid[t]) {
					// We can extend
					endTime = t;
				} else {
					terminateSequence = true;
				}
				if (t == var2Valid.length - 1) {
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
	 * @see infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate#getAddedMoreThanOneObservationSet()
	 */
	@Override
	public boolean getAddedMoreThanOneObservationSet() {
		return addedMoreThanOneObservationSet;
	}	

}
