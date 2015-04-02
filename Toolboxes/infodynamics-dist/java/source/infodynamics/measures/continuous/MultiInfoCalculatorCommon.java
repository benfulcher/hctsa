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

import java.util.Random;
import java.util.Vector;

/**
 * Implements {@link MultiInfoCalculator} to provide a base
 * class with common functionality for child class implementations of
 * {@link MultiInfoCalculator}
 * via various estimators. 
 * 
 * <p>These various estimators include: e.g. box-kernel estimation, KSG estimators, etc
 * (see the child classes linked above).
 * </p>
 * 
 * <p>Usage is as outlined in {@link MultiInfoCalculator}.</p>
 * 
  * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public abstract class MultiInfoCalculatorCommon implements MultiInfoCalculator {

	/**
	 * Number of joint variables to consider
	 */
	protected int dimensions = 0;
	/**
	 * Number of samples supplied
	 */
	protected int totalObservations = 0;
	/**
	 * Whether we are in debug mode
	 */
	protected boolean debug = false;
	/**
	 * Cached supplied observations
	 */
	protected double[][] observations;
	/**
	 * Set of individually supplied observations
	 */
	protected Vector<double[]> individualObservations;
	protected boolean miComputed = false;
	/**
	 * Cached last multi-info value calculated
	 */
	protected double lastAverage;
	/**
	 * Whether the user has supplied more than one (disjoint) set of samples
	 */
	protected boolean addedMoreThanOneObservationSet;

	/**
	 * Whether to normalise incoming values
	 */
	protected boolean normalise = true;

	private boolean underSample = false;
	private double samplingFactor = 0.1;
	private Random rand;

	@Override
	public void initialise(int dimensions) {
		this.dimensions = dimensions;
		lastAverage = 0.0;
		totalObservations = 0;
		miComputed = false;
		observations = null;
		addedMoreThanOneObservationSet = false;
	}

	/**
	 * Set properties for the calculator.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * <p>Valid property names, and what their
	 * values should represent, include properties defined by {@link MultiInfoCalculator#setProperty(String, String)}.</li>
	 * </ul>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
			normalise = Boolean.parseBoolean(propertyValue);
		} else if (propertyName.equalsIgnoreCase(SAMPLING_FACTOR_PROP_NAME)) {
			// Use less than 100 % of samples in the addObservation method
			samplingFactor = Double.parseDouble(propertyValue);
			underSample = (samplingFactor < 1.0);
			rand = new Random();
		} else {
			// No property was set here
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	@Override
	public void setObservations(double[][] observations) throws Exception {
		startAddObservations();
		addObservations(observations);
		finaliseAddObservations();
		addedMoreThanOneObservationSet = false;
	}

	@Override
	public void startAddObservations() {
		individualObservations = new Vector<double[]>();
	}

	@Override
	public void addObservation(double[] observation) {
		if (underSample && (rand.nextDouble() >= samplingFactor)) {
			// Don't take this sample
			return;
		}
		if (individualObservations.size() > 0) {
			addedMoreThanOneObservationSet = true;
		}
		individualObservations.add(observation);
	}

	@Override
	public void addObservations(double[][] observations) {
		// This implementation is not particularly efficient,
		//  however it will suffice for now.
		for (int s = 0; s < observations.length; s++) {
			addObservation(observations[s]);
		}
	}

	/**
	 * {@inheritDoc} 
	 * 
	 * This class provides a basic implementation, generating
	 * the internal set of samples in observations; child classes
	 * should then process these observations as required.
	 * 
	 */
	public void finaliseAddObservations() throws Exception {
		observations = new double[individualObservations.size()][];
		for (int t = 0; t < observations.length; t++) {
			observations[t] = individualObservations.elementAt(t);
		}
		// Allow vector to be reclaimed
		individualObservations = null;
		
		if (observations[0].length != dimensions) {
			throw new Exception("Incorrect number of dimensions " + observations[0].length +
					" in supplied observations (expected " + dimensions + ")");
		}
		totalObservations = observations.length;
	}
	
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		int[][][] newOrderings = new int[numPermutationsToCheck][][];
		// Generate numPermutationsToCheck * (dimensions-1) permutations of 0 .. data.length-1
		for (int n = 0; n < numPermutationsToCheck; n++) {
			// (Not necessary to check for distinct random perturbations)
			newOrderings[n] = rg.generateRandomPerturbations(totalObservations, dimensions-1);
		}
		return computeSignificance(newOrderings);
	}

	@Override
	public EmpiricalMeasurementDistribution computeSignificance(int[][][] newOrderings) throws Exception {
		
		int numPermutationsToCheck = newOrderings.length;
		if (!miComputed) {
			computeAverageLocalOfObservations();
		}
		
		// Store the real observations and their MI:
		double actualMI = lastAverage;
		
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		
		int countWhereMiIsMoreSignificantThanOriginal = 0;
		for (int i = 0; i < numPermutationsToCheck; i++) {
			// Compute the MI under this reordering
			double newMI = computeAverageLocalOfObservations(newOrderings[i]);
			measDistribution.distribution[i] = newMI;
			if (debug){
				System.out.println("New MI was " + newMI);
			}
			if (newMI >= actualMI) {
				countWhereMiIsMoreSignificantThanOriginal++;
			}
		}
		
		// Restore the actual MI and the observations
		lastAverage = actualMI;

		// And return the significance
		measDistribution.pValue = (double) countWhereMiIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualMI;
		return measDistribution;
	}

	@Override
	public double computeAverageLocalOfObservations(int[][] newOrdering)
			throws Exception {
		
		if (newOrdering == null) {
			return computeAverageLocalOfObservations();
		}
		
		// Take a clone of the object to compute the multi-info of the surrogates:
		// (this is a shallow copy, it doesn't make new copies of all
		//  the arrays)
		MultiInfoCalculatorCommon miSurrogateCalculator =
				(MultiInfoCalculatorCommon) this.clone();
		
		// Generate a new re-ordered source data
		double[][] shuffledData =
				MatrixUtils.reorderDataForVariables(
					observations, newOrdering);
		// Perform new initialisations
		miSurrogateCalculator.initialise(dimensions);
		// Set new observations
		miSurrogateCalculator.setObservations(shuffledData);
		// Compute the MI
		return miSurrogateCalculator.computeAverageLocalOfObservations();
	}

	@Override
	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	@Override
	public double getLastAverage() {
		return lastAverage;
	}
}
