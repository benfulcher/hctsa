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

package infodynamics.measures.mixed;

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.utils.EmpiricalMeasurementDistribution;

/**
 * <p>Interface defining computation of the differential mutual information between a given multivariate set of
 *  observations 
 *  and a discrete variable.
 *  This is done by examining the conditional probability distribution (given the discrete
 *  variable) against the probability distribution for the multivariate set.</p>
 *  
 * <p>
 * Usage intended for the child classes:
 * 	<ol>
 * 		<li>Construct the child class</li>
 *		<li>{@link #initialise(int, int)} to tell the class the number of dimensions of the continuous observations,
 *          and the base (number of available states) of the discrete variable.</li>
 * 		<li>Set properties using {@link #setProperty(String, String)}</li>
 * 		<li>Provide the observations to the calculator using:
 * 			{@link #setObservations(double[][], int[])}.</li>
 * 		<li>Compute the required information-theoretic results, primarily:
 * 			{@link #computeAverageLocalOfObservations()} to return the average 
 *          MI; or other calls to compute
 *          local values or statistical significance.</li>
 * 	</ol>
 * </p>
 * 
 * <p>These calculators are <b>EXPERIMENTAL</b> -- not properly tested,
 * and not well documented. The intended calling pattern is similar to
 * {@link MutualInfoCalculatorMultiVariate}
 * </p>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public interface MutualInfoCalculatorMultiVariateWithDiscrete {

	/**
	 * Initialise the calculator for use or reuse.
	 * This clears any previously supplied observations, but 
	 *  preserves the supplied properties.
	 * 
	 * @param dimensions number of joint continuous variables
	 * @param base number of states in the discrete observations
	 * @throws Exception
	 */
	public void initialise(int dimensions, int base) throws Exception;
	
	/**
	 * <p>Set the required property of the calculator to the given value.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * </p>
	 * 
	 * <p>There are no general properties settable on all child classes;
	 * each child class may define their own properties.</p>
	 * 
	 * @param propertyName name of property
	 * @param propertyValue value of property
	 */
	public void setProperty(String propertyName, String propertyValue);
	
	/**
	 * Set the properties from which the mutual information should be computed.
	 * 
	 * @param continuousObservations observations of the joint continuous variables;
	 *    first index is time or observation number, second index is variable number
	 *    (the number of variables should be equal to the number of dimensions set
	 *    in {@link #initialise(int, int)}.)
	 * @param discreteObservations observations of the discrete variable; must be
	 *    the same number of observations as supplied for continuousObservations.
	 *    Each observation must lie in the range 0..base-1, where base was set by
	 *    {@link #initialise(int, int)}.
	 * @throws Exception
	 */
	public void setObservations(double[][] continuousObservations,
			int[] discreteObservations) throws Exception;

	/**
	 * Compute the average mutual information from the previously supplied
	 *  observations via {@link #setObservations(double[][], int[])}
	 * 
	 * @return a scalar for the average mutual information
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations() throws Exception;

	/**
	 * Compute the local mutual information for each pair of continuous
	 *  and discrete values supplied here, using probability distribution
	 *  functions generated using the observations previously supplied
	 *  to the method {@link #setObservations(double[][], int[])}.
	 * 
	 * @param contStates observations of the joint continuous variables,
	 *  as specified in {@link #setObservations(double[][], int[])}
	 * @param discreteStates observations of the discrete variable,
	 *  as specified in {@link #setObservations(double[][], int[])}
	 * @return a time series of the local mutual information values
	 *  for each supplied observation
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] contStates, int[] discreteStates) throws Exception;

	/**
	 * <p>Compute the empiricial statistical significance of the mutual information
	 *  of the previously supplied observations
	 *  via {@link #setObservations(double[][], int[])}.
	 * We destroy the p(x,y) correlations, while retaining the p(x), p(y) marginals, to check how
	 *  significant this mutual information actually was.
	 * Specifically, this method returns an {@link EmpiricalMeasurementDistribution}
	 *  object containing numPermutationsToCheck surrogate measurements 
	 *  where the discrete data is shuffled against the continuous data set.
	 * </p>
	 *  
	 * <p>This is in the spirit of Chavez et. al.
	 *  which was performed for Transfer entropy.
	 * </p>
	 * 
	 * @param numPermutationsToCheck the number of permuted surrogates to examine
	 * @return the proportion of MI scores from the distribution which have higher or equal MIs to ours.
	 * @see "Chavez et. al., 'Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals', Journal of Neuroscience Methods 124 (2003) 113-128"
	 * @throws Exception
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception;

	/**
	 * <p>Compute the empiricial statistical significance of the mutual information
	 *  of the previously supplied observations
	 *  via {@link #setObservations(double[][], int[])}.
	 *  This method performs as per {@link #computeSignificance(int)} except
	 *  that the shuffling of the discrete time series is 
	 *  specified here by the newOrderings parameter.
	 * </p>
	 * 
	 * @param newOrderings the specific new orderings to use. First index
	 *  is for the new ordering number; second index is, for that
	 *  particular reordering, which time point of the original series
	 *  to grab at that point. 
	 * @return the proportion of MI scores from the distribution which have higher or equal MIs to ours.
	 * @see #computeSignificance(int)
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) throws Exception;

	/**
	 * Set whether to print extra debug messages
	 * 
	 * @param debug whether to print extra debug messages
	 */
	public void setDebug(boolean debug);

	/**
	 * Return the previously computed average
	 * 
	 * @return the previously computed average
	 */
	public double getLastAverage();
	
	/**
	 * Get the number of observations supplied via
	 *  {@link #setObservations(double[][], int[])}
	 * 
	 * @return the number of supplied observations
	 */
	public int getNumObservations();

}
