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

package infodynamics.measures.mixed.kernel;

import java.util.Arrays;

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.kernel.KernelEstimatorMultiVariate;
import infodynamics.measures.mixed.MutualInfoCalculatorMultiVariateWithDiscrete;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

/**
 * <p>These calculators are <b>EXPERIMENTAL</b> -- not properly tested,
 * and not well documented. The intended calling pattern is similar to
 * {@link MutualInfoCalculatorMultiVariate}
 * </p>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class MutualInfoCalculatorMultiVariateWithDiscreteKernel implements
	MutualInfoCalculatorMultiVariateWithDiscrete {

	protected KernelEstimatorMultiVariate mvke = null;
	protected KernelEstimatorMultiVariate[] mvkeForEachDiscrete = null;
	protected int base = 0;

	private int totalObservations = 0;
	// private int dimensions1 = 0;
	// private int dimensions2 = 0;
	private boolean debug = false;
	private double[][] contObservations; 
	private int[] discObservations;
	private int[] discCounts;
	private double lastAverage;
	private boolean miComputed;
	
	private boolean normalise = true;
	public static final String NORMALISE_PROP_NAME = "NORMALISE";

	// No dynamic correlation exclusion time, since we won't track time
	//  for the conditional observations
	
	private boolean forceCompareToAll = false;
	public static final String FORCE_KERNEL_COMPARE_TO_ALL = "FORCE_KERNEL_COMPARE_TO_ALL";
	
	/**
	 * Default value for epsilon
	 */
	public static final double DEFAULT_EPSILON = 0.25;
	/**
	 * Kernel width
	 */
	private double epsilon = DEFAULT_EPSILON;
	private boolean usingSingleKernelWidthValue = true;
	private double[] epsilons = null;
	public static final String EPSILON_PROP_NAME = "EPSILON";
	
	public MutualInfoCalculatorMultiVariateWithDiscreteKernel() {
		mvke = new KernelEstimatorMultiVariate();
		mvke.setNormalise(normalise);
	}

	/**
	 * Initialise using a the current settings for the kernel width
	 *  (which is the default kernel width {@link DEFAULT_EPSILON}
	 *   if it has not yet been set)
	 * 
	 * @param dimensions for number of continuous variables
	 * @param base for discrete variable
	 */
	public void initialise(int dimensions, int base) {
		if (usingSingleKernelWidthValue) {
			mvke.initialise(dimensions, epsilon);
		} else {
			// epsilons must have been initialised previously in this case
			mvke.initialise(epsilons);
		}
		initialiseCommon(base);
	}

	/**
	 * Initialise using the supplied kernel width for all continuous variables
	 * 
	 * @param dimensions for number of continuous variables
	 * @param base for discrete variable
	 * @param epsilon kernel width
	 */
	public void initialise(int dimensions, int base, double epsilon) {
		this.epsilon = epsilon;
		usingSingleKernelWidthValue = true;
		mvke.initialise(dimensions, epsilon);
		initialiseCommon(base);
	}

	/**
	 * Initialise using the supplied kernel width for all continuous variables
	 * 
	 * @param base for discrete variable
	 * @param epsilons kernel width for each continuous variable
	 */
	public void initialise(int base, double epsilons[]) {
		this.epsilons = epsilons;
		usingSingleKernelWidthValue = false;
		mvke.initialise(epsilons);
		initialiseCommon(base);
	}

	protected void initialiseCommon(int base) {
		this.base = base;
		mvkeForEachDiscrete = new KernelEstimatorMultiVariate[base];
		for (int i = 0; i < base; i++) {
			mvkeForEachDiscrete[i] = new KernelEstimatorMultiVariate();
			// We won't normalise the conditional calculators, but feed them
			//  the same kernel width used by the full marginal one
			mvkeForEachDiscrete[i].setNormalise(false);
			mvkeForEachDiscrete[i].setForceCompareToAll(forceCompareToAll);
			// Don't initialise these calculators yet - we'll wait until
			//  we have the adjusted epsilon from the full marginal space
		}
		discCounts = new int[base];
		// this.dimensions1 = dimensions1;
		// this.dimensions2 = dimensions2;
		lastAverage = 0.0;
		miComputed = false;
	}
	
	/**
	 * Set the observations for the PDFs.
	 * Should only be called once, the last call contains the
	 *  observations that are used (they are not accumulated). 
	 * 
	 * @param observations
	 */
	public void setObservations(double continuousObservations[][], int discreteObservations[]) throws Exception {
		if (continuousObservations.length != discreteObservations.length) {
			throw new Exception("Observations are not of the same length");
		}
		this.contObservations = continuousObservations;
		mvke.setObservations(continuousObservations);
		setDiscreteData(continuousObservations, discreteObservations);
		totalObservations = continuousObservations.length;
	}
	
	protected void setDiscreteData(double continuousObservations[][], int discreteObservations[]) {
		// Clear the discrete counts:
		Arrays.fill(discCounts, 0);
		// Compute the observation counts for the discrete state
		for (int t = 0; t < discreteObservations.length; t++) {
			discCounts[discreteObservations[t]]++;
		}
		for (int i = 0; i < base; i++) {
			// Extract the observations for when this base value occurs:
			double[][] obsForThisDiscValue = MatrixUtils.extractSelectedPointsMatchingCondition(
					continuousObservations, discreteObservations, i, discCounts[i]);
			// Set the kernel width for the relevant kernel estimator:
			mvkeForEachDiscrete[i].initialise(mvke.getKernelWidthsInUse());
			// Set these observations for the relevant kernel estimator:
			mvkeForEachDiscrete[i].setObservations(obsForThisDiscValue);
		}
		this.discObservations = discreteObservations;
	}
	
	/**
	 * Compute the MI from the observations we were given.
	 * 
	 * @return MI in bits
	 */
	public double computeAverageLocalOfObservations() {
		double mi = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double probCont = mvke.getProbability(contObservations[b]);
			double condProbCont = mvkeForEachDiscrete[discObservations[b]].getProbability(contObservations[b]);
			double logTerm = 0.0;
			double cont = 0.0;
			if (condProbCont > 0.0) {
				// If we have counted joint correlations, we must have marginals for each
				logTerm = condProbCont / probCont;
				cont = Math.log(logTerm);
			}
			mi += cont;
			if (debug) {
				System.out.printf("%d: %.3f, %d, (%.3f %d, %.5f %d) %.5f -> %.5f -> %.5f\n",
						b, contObservations[b][0], discObservations[b],
						condProbCont, mvkeForEachDiscrete[discObservations[b]].getCount(contObservations[b]),
						probCont, mvke.getCount(contObservations[b]),
						logTerm, cont, mi);
			}
		}
		lastAverage = mi / (double) totalObservations / Math.log(2.0);
		miComputed = true;
		return lastAverage;
	}

	/**
	 * Compute the MI if data were reordered.
	 * 
	 * @param newOrdering
	 * @return MI under the reordering scheme
	 */
	public double computeAverageLocalOfObservations(int[] newOrdering) throws Exception {
		// Store the real observations and their MI:
		double actualMI = lastAverage;
		int[] originalDiscrete = discObservations;
		
		// Generate a new re-ordered data2
		int[] newDiscrete = MatrixUtils.extractSelectedTimePoints(originalDiscrete, newOrdering);
		
		// Perform new initialisations on the discrete pdfs
		setDiscreteData(contObservations, newDiscrete);
		// Compute the MI
		double newMI = computeAverageLocalOfObservations();
		
		// Restore the actual MI and the observations
		lastAverage = actualMI;
		setDiscreteData(contObservations, originalDiscrete);

		return newMI;
	}

	public synchronized EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(contObservations.length, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}

	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) throws Exception {
		int numPermutationsToCheck = newOrderings.length;
		if (!miComputed) {
			computeAverageLocalOfObservations();
		}
		// Store the real observations and their MI:
		double actualMI = lastAverage;
		int[] originalDiscrete = discObservations;
		
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		
		int countWhereMiIsMoreSignificantThanOriginal = 0;
		for (int i = 0; i < numPermutationsToCheck; i++) {
			// Generate a new re-ordered data2
			int[] newDiscrete = MatrixUtils.extractSelectedTimePoints(originalDiscrete, newOrderings[i]);
			// Perform new initialisations on the discrete pdfs
			setDiscreteData(contObservations, newDiscrete);
			// Compute the MI
			double newMI = computeAverageLocalOfObservations();
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
		setDiscreteData(contObservations, originalDiscrete);

		// And return the significance
		measDistribution.pValue = (double) countWhereMiIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualMI;
		
		return measDistribution;
	}

	/**
	 * Extra utility method to return the joint entropy
	 * 
	 * @return
	 */
	public double computeAverageJointEntropy() {
		double entropy = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob = mvkeForEachDiscrete[discObservations[b]].getProbability(contObservations[b])
							* (double) discCounts[discObservations[b]] / (double) totalObservations;
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
	 * Extra utility method to return the entropy of the first set of joint variables
	 * 
	 * @return
	 */
	public double computeAverageEntropyOfObservation1() {
		double entropy = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob = mvke.getProbability(contObservations[b]);
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
	 * Extra utility method to return the entropy of the second set of joint variables
	 * 
	 * @return
	 */
	public double computeAverageEntropyOfObservation2() {
		double entropy = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob = (double) discCounts[discObservations[b]] / (double) totalObservations;
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
	 * Extra utility method to return the information distance
	 * 
	 * @return
	 */
	public double computeAverageInfoDistanceOfObservations() {
		throw new RuntimeException("Not implemented yet");
		/*
		double infoDistance = 0.0;
		for (int b = 0; b < totalObservations; b++) {
			double prob1 = mvke.getProbability(observations1[b], b);
			double prob2 = mvke2.getProbability(observations2[b], b);
			double probJoint = mvkeJoint.getProbability(observations1[b], observations2[b], b);
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
		*/
	}

	/**
	 * Compute the local MI values for the previous observations.
	 * 
	 * @return
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		return computeLocalUsingPreviousObservations(contObservations, discObservations);
	}

	/**
	 * Compute the local MI values for these given values, using the previously provided
	 *  observations to compute the probabilities.
	 * 
	 * @param states1
	 * @param states2
	 * @return
	 */
	public double[] computeLocalUsingPreviousObservations(double states1[][],
			int[] states2) {
		
		double mi = 0.0;
		int timeSteps = states1.length;
		double[] localMi = new double[timeSteps];
		double condProbCont, probCont;
		for (int b = 0; b < timeSteps; b++) {
			probCont = mvke.getProbability(states1[b]);
			condProbCont = mvkeForEachDiscrete[states2[b]].getProbability(states1[b]);
			double logTerm = 0.0;
			localMi[b] = 0.0;
			if (condProbCont > 0.0) {
				// By necessity prob1 and prob2 will be > 0.0
				logTerm = condProbCont / probCont;
				localMi[b] = Math.log(logTerm) / Math.log(2.0);
			}
			mi += localMi[b];
			if (debug) {
				System.out.printf("%d: (%.5f, %.5f) %.5f -> %.5f -> %.5f\n",
						b, condProbCont, probCont, logTerm, localMi[b], mi);
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
	 * @param states1
	 * @param states2
	 * @return
	 */
	public double[] computeLocalJointEntropyOfPreviousObservations() throws Exception {
		return computeLocalJointEntropyUsingPreviousObservations(contObservations,
				discObservations);
	}

	/**
	/**
	 * Internal implementation
	 * 
	 * @param states1
	 * @param states2
	 * @param isOurPreviousObservations
	 * @return
	 */
	public double[] computeLocalJointEntropyUsingPreviousObservations(
			double states1[][], int states2[]) {
		int timeSteps = states1.length;
		double[] localJoint = new double[timeSteps];
		double prob;
		for (int b = 0; b < totalObservations; b++) {
			prob = mvkeForEachDiscrete[states2[b]].getProbability(states1[b]) *
				(double) discCounts[states2[b]] / (double) totalObservations;
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
	 *  observations for VARIABLE 1 to compute the probabilities.
	 * 
	 * @param states1
	 *  
	 * @return
	 */
	public double[] computeLocalEntropy1OfPreviousObservations() {
		return computeLocalEntropyFromPreviousObservations(contObservations);
	}


	/**
	 * Compute the local entropy values for the previously provided
	 *  observations for VARIABLE 2 to compute the probabilities.
	 * 
	 * @param states2
	 *  
	 * @return
	 */
	public double[] computeLocalEntropy2OfPreviousObservations() {
		return computeLocalEntropyFromPreviousObservations(discObservations);
	}

	/**
	 * Utility function to implement computeLocalEntropy1FromPreviousObservations
	 *  
	 * @param states
	 * @return
	 */
	public double[] computeLocalEntropyFromPreviousObservations(
			double states[][]) {
		int timeSteps = states.length;
		double[] localEntropy = new double[timeSteps];
		double prob;
		for (int b = 0; b < totalObservations; b++) {
			prob = mvke.getProbability(states[b]);
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
	 * Utility function to implement computeLocalEntropy1FromPreviousObservations
	 *  
	 * @param states
	 * @return
	 */
	public double[] computeLocalEntropyFromPreviousObservations(
			int states[]) {
		int timeSteps = states.length;
		double[] localEntropy = new double[timeSteps];
		double prob;
		for (int b = 0; b < totalObservations; b++) {
			prob = (double) discCounts[states[b]] / (double) totalObservations;
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

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public double getLastAverage() {
		return lastAverage;
	}

	/**
	 * <p>Set properties for the mutual information calculator.
	 * These can include:
	 * <ul>
	 * 		<li>{@link #EPSILON_PROP_NAME} - applies to full marginal space of continuous</li>
	 * 		<li>{@link #NORMALISE_PROP_NAME}</li>
	 * 		<li>{@link #FORCE_KERNEL_COMPARE_TO_ALL}</li>
	 * </ul></p>
	 * 
	 * <p>Note that dynamic correlation exclusion may have unexpected results if multiple
	 *  observation sets have been added. This is because multiple observation sets
	 *  are treated as though they are from a single time series, so observations from
	 *  near the end of observation set i will be excluded from comparison to 
	 *  observations near the beginning of observation set (i+1). 
	 * </p>
	 * 
	 * @param propertyName
	 * @param propertyValue
	 */
	public void setProperty(String propertyName, String propertyValue) {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(EPSILON_PROP_NAME)) {
			usingSingleKernelWidthValue = true;
			epsilon = Double.parseDouble(propertyValue);
		} else if (propertyName.equalsIgnoreCase(NORMALISE_PROP_NAME)) {
			normalise = Boolean.parseBoolean(propertyValue);
			mvke.setNormalise(normalise);
			// Don't set the normalise property on the conditional
			//  kernel estimation calculators - we'll set these directly
			//  from the joint space
		} else if (propertyName.equalsIgnoreCase(FORCE_KERNEL_COMPARE_TO_ALL)) {
			forceCompareToAll = Boolean.parseBoolean(propertyValue);
			mvke.setForceCompareToAll(forceCompareToAll);
			for (int i = 0; i < base; i++) {
				mvkeForEachDiscrete[i].setForceCompareToAll(forceCompareToAll);
			}
		} else {
			// No property was set
			propertySet = false;
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	public int getNumObservations() {
		return totalObservations;
	}
	
	/**
	 * Return the kernel widths used by the MI calculator for the
	 *  continuous variables here.
	 * This should not be called until the observations have been set.
	 * These are the kernel widths actually applied to the data (not the 
	 *  number of standard deviations, if we are using normalisation).
	 * 
	 * @return an array of doubles with the kernel widths.
	 */
	public double[] getKernelWidthsInUse() {
		// Return a copy so that the user can't mess with it
		return Arrays.copyOf(mvke.getKernelWidthsInUse(), mvke.getKernelWidthsInUse().length);
	}
}
