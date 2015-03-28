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

package infodynamics.measures.mixed.gaussian;

import infodynamics.measures.mixed.ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSourceCommon;
import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateCommon;
import infodynamics.measures.continuous.gaussian.EntropyCalculatorMultiVariateGaussian;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential conditional mutual information of a given multivariate set of
 *  observations with a discrete variable, conditioned on another multivariate set
 *  of observations,
 *  assuming that the probability distribution function for these continuous observations is
 *  a multivariate Gaussian distribution.</p>
 *  
 * <p>This is done by examining the conditional probability distribution for
 *  the multivariate continuous variable C (given the discrete
 *  variable D) against the probability distribution for C, all conditioned
 *  on another continuous variable Z:
 *  MI(C;D|Z) := H(C|Z) - H(C|D,Z).</p>
 *  
 * <p><b>CAVEAT EMPTOR</b>: The real question this type of calculation asks
 * is to what extent does knowing the value of the discrete variable reduce
 * variance in the continuous variable(s), given the other known continuous variable.
 *  Indeed, it may not to behave (I should test this--)
 * as we would normally expect a mutual information calculation: if we add more
 * continuous variables in, it may increase in spite of redundancy between these
 * variables. TODO Further exploration should take place here ...</p>
 * 
 * <p>These calculators are <b>EXPERIMENTAL</b> -- not properly tested,
 * and not well documented. The intended calling pattern is similar to
 * {@link ConditionalMutualInfoCalculatorMultiVariate}
 * </p>
 * 
 * <p>
 * Usage:
 * 	<ol>
 * 		<li>Construct {@link #ConditionalMutualnfoCalculatorMultiVariateWithDiscreteSourceGaussian()}</li>
 *		<li>{@link #initialise(int, int, int)}</li>
 * 		<li>Set properties using {@link #setProperty(String, String)}</li>
 * 		<li>Provide the observations to the calculator using:
 * 			{@link #setObservations(double[][], int[], double[][])}, or
 * 			a sequence of:
 * 			{@link #startAddObservations()},
 *          multiple calls to {@link #addObservations(double[][], int[], double[][])}
 *          and then
 *          {@link #finaliseAddObservations()}.</li>
 * 		<li>Compute the required information-theoretic results, primarily:
 * 			{@link #computeAverageLocalOfObservations()} to return the average differential
 *          entropy based on the variance of
 *          the supplied observations; or other calls to compute
 *          local values or statistical significance.</li>
 * 	</ol>
 * </p>
 * 
 * <p>
 * Alters behaviour slightly from parent class {@link ConditionalMutualInfoMultiVariateCommon}
 * in that property {@link ConditionalMutualInfoMultiVariateCommon#PROP_NORMALISE}
 * is set to false by default here (since this makes more sense for
 * linear-Gaussian analysis). 
 * </p>
 * 
 * @see <a href="http://mathworld.wolfram.com/DifferentialEntropy.html">Differential entropy for Gaussian random variables at Mathworld</a>
 * @see <a href="http://en.wikipedia.org/wiki/Differential_entropy">Differential entropy for Gaussian random variables at Wikipedia</a>
 * @see <a href="http://en.wikipedia.org/wiki/Multivariate_normal_distribution">Multivariate normal distribution on Wikipedia</a>
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSourceGaussian
		extends	ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSourceCommon 
		implements Cloneable {

	/**
	 * Entropy calculator applied to the whole set of conditional data Z
	 */
	protected EntropyCalculatorMultiVariateGaussian entCalcZ;
	/**
	 * Entropy calculator applied to the whole set of continuous data C
	 * and conditional data Z
	 * 
	 */
	protected EntropyCalculatorMultiVariateGaussian entCalcCZ;
	/**
	 * Entropy calculators applied to the set of conditional data Z
	 * associated with each discrete value
	 */
	protected EntropyCalculatorMultiVariateGaussian[] entCalcZForEachDiscrete;
	/**
	 * Entropy calculators applied to the set of continuous data C
	 * and conditional data Z
	 * associated with each discrete value
	 */
	protected EntropyCalculatorMultiVariateGaussian[] entCalcCZForEachDiscrete;

	public ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSourceGaussian() {
		super();
		// Normalising data makes less sense for linear-Gaussian estimation,
		//  so we turn this off by default.
		normalise = false;
		entCalcZ = new EntropyCalculatorMultiVariateGaussian();
		entCalcZForEachDiscrete = null;
		entCalcCZ = new EntropyCalculatorMultiVariateGaussian();
		entCalcCZForEachDiscrete = null;
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSourceCommon#initialise(int, int, int)
	 */
	@Override
	public void initialise(int dimensions, int base, int dimensionsCond) {
		super.initialise(dimensions, base, dimensionsCond);
		entCalcZ.initialise(dimensionsCond);
		entCalcCZ.initialise(dimensions + dimensionsCond);
		entCalcZForEachDiscrete = new EntropyCalculatorMultiVariateGaussian[base];
		entCalcCZForEachDiscrete = new EntropyCalculatorMultiVariateGaussian[base];
		for (int b = 0; b < base; b++) {
			entCalcZForEachDiscrete[b] = new EntropyCalculatorMultiVariateGaussian();
			entCalcCZForEachDiscrete[b] = new EntropyCalculatorMultiVariateGaussian();
			// If any properties relevant for these calculators were set in
			//  setProperty then we should set them here
			entCalcZForEachDiscrete[b].initialise(dimensionsCond);
			entCalcCZForEachDiscrete[b].initialise(dimensions + dimensionsCond);
		}
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSourceCommon#finaliseAddObservations()
	 */
	@Override
	public void finaliseAddObservations() throws Exception {
		super.finaliseAddObservations();
		
		// Set the complete set of observations:
		//  (this will pick up any errors in the dimensions of the continuous
		//   observations)
		entCalcZ.setObservations(conditionedDataZ);
		double[][] joinedCZ = MatrixUtils.appendColumns(continuousDataX, conditionedDataZ);
		entCalcCZ.setObservations(joinedCZ);
		
		// Set the observations corresponding to each discrete value:
		int totalNumberOfSuppliedObservations = 0;
		for (int b = 0; b < base; b++) {
			// Extract the observations for when this base value occurs:
			double[][] obsZForThisDiscValue = MatrixUtils.extractSelectedPointsMatchingCondition(
					conditionedDataZ, discreteData, b);
			double[][] obsCZForThisDiscValue = MatrixUtils.extractSelectedPointsMatchingCondition(
					joinedCZ, discreteData, b);
			// Set the observations for each discrete value:
			entCalcZForEachDiscrete[b].setObservations(obsZForThisDiscValue);
			entCalcCZForEachDiscrete[b].setObservations(obsCZForThisDiscValue);
			totalNumberOfSuppliedObservations += obsZForThisDiscValue.length;
		}
		// Check that all of the supplied observations were extracted corresponding
		//  to one of the allowed discrete values
		if (totalNumberOfSuppliedObservations != discreteData.length) {
			throw new Exception("Some values in discreteObservations were not in the range 0..base-1");
		}
	}

	public double computeAverageLocalOfObservations() throws Exception {
		// The average mutual information can be expressed
		//  as a difference between the conditional entropy of the continuous observations
		//  and the conditional entropy of the continuous given the
		//  discrete observations, both given the conditional continuous observations:
		//		I(C;D|Z) = H(C|Z) - H(C|D,Z)
		//               = H(C,Z) - H(Z) - H(C,Z|D) + H(C|D)
		double meanConditionalEntropyCZ = 0;
		double meanConditionalEntropyZ = 0;
		for (int b = 0; b < base; b++) {
			double pOfB = (double) counts[b] / (double) totalObservations;
			meanConditionalEntropyZ += pOfB *
					entCalcZForEachDiscrete[b].computeAverageLocalOfObservations();
			meanConditionalEntropyCZ += pOfB *
					entCalcCZForEachDiscrete[b].computeAverageLocalOfObservations();
		}
		double entCZ = entCalcCZ.computeAverageLocalOfObservations();
		double entZ = entCalcZ.computeAverageLocalOfObservations();
		condMi = entCZ - entZ
				- meanConditionalEntropyCZ
				+ meanConditionalEntropyZ;
		if (debug) {
			System.out.printf("H(C,Z)=%.4f - H(Z)=%.4f - H(C,Z|D)=%.4f + H(Z|D)=%.4f = %.4f\n",
					entCZ, entZ, meanConditionalEntropyCZ, meanConditionalEntropyZ, condMi);
		}
		return condMi;
	}

	public double[] computeLocalOfPreviousObservations() throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public double[] computeLocalUsingPreviousObservations(
			double[][] contStates, int[] discreteStates,
			double[][] conditionedStates) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	/**
	 * Clone the object - note: while it does create new cloned instances of
	 *  the {@link EntropyCalculatorMultiVariateGaussian} objects,
	 *  I think these only
	 *  have shallow copies to the data.
	 * This is enough though to maintain the structure across 
	 *  various {@link #computeSignificance(boolean, int)} calls.
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	protected Object clone() throws CloneNotSupportedException {
		ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSourceGaussian theClone = 
			(ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSourceGaussian) super.clone();
		// Now assign clones of the EntropyCalculatorMultiVariateGaussian objects:
		// First clone those for the conditional variable:
		theClone.entCalcZ =
				(EntropyCalculatorMultiVariateGaussian) entCalcZ.clone();
		if (entCalcZForEachDiscrete != null) {
			theClone.entCalcZForEachDiscrete = new EntropyCalculatorMultiVariateGaussian[base];
			for (int b = 0; b < base; b++) {
				theClone.entCalcZForEachDiscrete[b] =
						(EntropyCalculatorMultiVariateGaussian)
							entCalcZForEachDiscrete[b].clone();
			}
		}
		// Next clone those for the conditional variable and the continuous source variable
		theClone.entCalcCZ =
				(EntropyCalculatorMultiVariateGaussian) entCalcCZ.clone();
		if (entCalcCZForEachDiscrete != null) {
			theClone.entCalcCZForEachDiscrete = new EntropyCalculatorMultiVariateGaussian[base];
			for (int b = 0; b < base; b++) {
				theClone.entCalcCZForEachDiscrete[b] =
						(EntropyCalculatorMultiVariateGaussian)
							entCalcCZForEachDiscrete[b].clone();
			}
		}
		return theClone;
	}
}
