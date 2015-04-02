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

/**
 * Interface for implementations of
 * Active Information Storage estimators on continuous univariate data (ie double[] arrays).
 * 
 * <p>See definition of Active Information Storage (AIS) by Lizier et al. below.
 * Basically, AIS is the mutual information between the past <i>state</i>
 * of a time-series process <i>X</i> and its next value. The past <i>state</i> at time <code>n</code>
 * is represented by an embedding vector of <code>k</code> values from <code>X_n</code> backwards,
 * each separated by <code>\tau</code> steps, giving
 * <code><b>X^k_n</b> = [ X_{n-(k-1)\tau}, ... , X_{n-\tau}, X_n]</code>.
 * We call <code>k</code> the embedding dimension, and <code>\tau</code>
 * the embedding delay.
 * AIS is then the mutual information between <b>X^k_n</b> and X_{n+1}.</p>
 * 
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * 	<ol>
 * 		<li>Construct the calculator;</li>
 * 		<li>Set properties using {@link #setProperty(String, String)};</li>
 *		<li>Initialise the calculator using {@link #initialise()} or
 *			{@link #initialise(int)} or {@link #initialise(int, int)};
 *		</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using:
 * 			<ul>
 * 				<li>{@link #setObservations(double[])} or
 * 					{@link #setObservations(double[], boolean[])} for calculations
 * 					based on single time-series, OR</li>
 * 				<li>The following sequence:<ol>
 * 						<li>{@link #startAddObservations()}, then</li>
 * 						<li>One or more calls to {@link #addObservations(double[])} or
 * 							{@link #addObservations(double[], int, int)}, then</li>
 * 						<li>{@link #finaliseAddObservations()};</li>
 * 					</ol></li>
 * 			</ul>
 * 		</li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average AIS: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local AIS values for these samples: {@link #computeLocalOfPreviousObservations()}</li>
 * 				<li>local AIS values for a specific set of samples: {@link #computeLocalUsingPreviousObservations(double[])}</li>
 * 				<li>the distribution of AIS values under the null hypothesis
 * 					of no relationship between past sequences in the series
 * 					and the next value: {@link #computeSignificance(int)} or
 * 					{@link #computeSignificance(int[][])}.</li>
 * 			</ul>
 * 		</li>
 * 		<li>
 * 		Return to step 2 or 3 to re-use the calculator on a new data set.
 * 		</li>
 * 	</ol>
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
 */
public interface ActiveInfoStorageCalculator {

	/**
	 * Property name for embedding length <code>k</code> of
	 * the past history vector (1 by default).
	 */
	public static final String K_PROP_NAME = "k_HISTORY";
	/**
	 * Property name for embedding delay <code>\tau</code> of the past history vector
	 * (1 by default).
	 * The delay exists between points in the past k-vector
	 *  but there is always only one time step between the
	 *  past k-vector and the next observation.  
	 */
	public static final String TAU_PROP_NAME = "TAU";

	/**
	 * Initialise the calculator for (re-)use, with the existing (or default) values of parameters
	 * Clears any PDFs of previously supplied observations.
	 * 
	 */
	public void initialise() throws Exception;

	/**
	 * Initialise the calculator for (re-)use, with some parameters
	 * supplied here, and existing (or default) values of other parameters
	 * to be used.
	 * 
	 * @param k embedding length of past history vector
	 */
	public void initialise(int k) throws Exception;
	
	/**
	 * Initialise the calculator for (re-)use, with some parameters
	 * supplied here, and existing (or default) values of other parameters
	 * to be used.
	 * 
	 * @param k embedding length of past history vector
	 * @param tau embedding delay of past history vector
	 */
	public void initialise(int k, int tau) throws Exception;

	/**
	 * Set properties for the underlying calculator implementation.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * <p>Property names defined at the interface level, and what their
	 * values should represent, include:</p>
	 * <ul>
	 * 	<li>{@link #K_PROP_NAME} -- embedding length <code>k</code> of
	 * 		the past history vector</li>
	 *  <li>{@link #TAU_PROP_NAME} -- embedding delay between each of the
	 *     <code>k</code> points in the past history vector.</li>
	 * </ul>
	 *  
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * <p>Note that implementing classes may defined additional properties.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception;

	/**
	 * Sets a single time-series from which to compute the PDF for the AIS.
	 * Cannot be called in conjunction with other methods for setting/adding
	 * observations.
	 * 
	 * @param observations time-series array of (univariate) samples,
	 *  where the array index is time.
	 */
	public void setObservations(double observations[]) throws Exception;
	
	/**
	 * Signal that we will add in the samples for computing the PDF 
	 * from several disjoint time-series or trials via calls to
	 * {@link #addObservations(double[])} etc.
	 *
	 */
	public void startAddObservations();
	
	/**
	 * Add more time-series for the computation of the PDF.
	 * The array observations must not be over-written by the user
	 *  until after finaliseAddObservations() has been called.
	 * 
	 * @param observations time-series array of (univariate) samples,
	 *  where the array index is time.
	 */
	public void addObservations(double[] observations) throws Exception;

	/**
	 * Add more time-series for the computation of the PDF, using
	 * only a sub-series of <code>observations</code>.
	 * The array observations must not be over-written by the user
	 *  until after finaliseAddObservations() has been called.
	 * 
	 * @param observations time-series array of (univariate) samples,
	 *  where the array index is time.
	 * @param startTime first time index to extract samples from
	 * @param numTimeSteps number of time steps to extract starting from startTime
	 */
	public void addObservations(double[] observations,
			int startTime, int numTimeSteps) throws Exception ;

	/**
	 * Signal that the observations are now all added, PDFs can now be constructed.
	 * 
	 * @throws Exception 
	 */
	public void finaliseAddObservations() throws Exception;
	
	/**
	 * Sets a single time-series from which to compute the PDF for the AIS,
	 * subject to the validity of each sample in that series.
	 * Cannot be called in conjunction with other methods for setting/adding
	 * observations.
	 *  
	 * @param observations time-series array of (univariate) samples,
	 *  where the array index is time.
	 * @param valid a time series (with indices the same as observations)
	 *  indicating whether the entry in observations at that index is valid; 
	 *  we only take vectors as samples to add to the observation set where
	 *  all points in the time series (even between points in 
	 *  the embedded k-vector with embedding delays) are valid.
	 */
	public void setObservations(double[] observations,
			boolean[] valid) throws Exception;

	/**
	 * Compute the AIS from the previously-supplied samples.
	 * 
	 * @return the AIS estimate
	 */
	public double computeAverageLocalOfObservations() throws Exception;

	/**
	 * Compute the local AIS values for each of the
	 * previously-supplied samples.
	 * 
	 * <p>PDFs are computed using all of the previously supplied
	 * observations.</p>
	 * 
	 * <p>If the samples were supplied via a single call such as
	 * {@link #setObservations(double[])},
	 * then the return value is a single time-series of local
	 * AIS values corresponding to these samples (with the first
	 * <code>(k-1)*tau + 1</code>
	 * values set to 0 since AIS is undefined there).</p>
	 * 
	 * <p>Otherwise where disjoint time-series observations were supplied using several 
	 *  calls such as {@link addObservations(double[])}
	 *  then the local values for each disjoint observation set will be appended here
	 *  to create a single "time-series" return array (without any <code>(k-1)*tau + 1</code>
	 *  leading 0 values).</p>
	 *  
	 * @return the "time-series" of local AIS values.
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception;

	/**
	 * Compute the local AIS values for each of the
	 * supplied samples in <code>newObservations</code>.
	 * 
	 * <p>PDFs are computed using all of the previously supplied
	 * observations, but not those in <code>newObservations</code> (unless they were
	 * some of the previously supplied samples).</p>
	 * 
	 * @param newObservations time-series for which to compute local AIS values
	 * @return time-series of local AIS values corresponding to newObservations
	 *   (the first <code>(k-1)*tau + 1</code>
	 *    values are set to 0 since AIS is undefined there)
	 * @throws Exception for an invalid input array, e.g. of length less than <code>k+1</code>
	 */
	public double[] computeLocalUsingPreviousObservations(double[] newObservations) throws Exception;

	/**
	 * Generate a bootstrapped distribution of what the AIS would look like,
	 * under a null hypothesis that the previous <code>k</code> values of our
	 * samples had no relation to the next value in the time-series.
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for AIS 
	 * as a mutual information. Basically, the marginal PDFs
	 * of the past <code>k</code> values, and that of the next value, 
	 * are preserved, while their joint PDF is destroyed, and the 
	 * distribution of AIS under these conditions is generated.</p>
	 * 
	 * <p>Note that if several disjoint time-series have been added 
	 * as observations using {@link #addObservations(double[])} etc.,
	 * then these separate "trials" will be mixed up in the generation
	 * of surrogates here.</p>
	 * 
	 * <p>This method (in contrast to {@link #computeSignificance(int[][])})
	 * creates <i>random</i> shufflings of the next values for the surrogate AIS
	 * calculations.</p>
	 * 
	 * @param numPermutationsToCheck number of surrogate samples to bootstrap
	 *  to generate the distribution.
	 * @return the distribution of AIS scores under this null hypothesis.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception;
	
	/**
	 * Generate a bootstrapped distribution of what the AIS would look like,
	 * under a null hypothesis that the previous <code>k</code> values of our
	 * samples had no relation to the next value in the time-series.
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for AIS 
	 * as a mutual information. Basically, the marginal PDFs
	 * of the past <code>k</code> values, and that of the next value, 
	 * are preserved, while their joint PDF is destroyed, and the 
	 * distribution of AIS under these conditions is generated.</p>
	 * 
	 * <p>Note that if several disjoint time-series have been added 
	 * as observations using {@link #addObservations(double[])} etc.,
	 * then these separate "trials" will be mixed up in the generation
	 * of surrogates here.</p>
	 * 
	 * <p>This method (in contrast to {@link #computeSignificance(int)})
	 * allows the user to specify how to construct the surrogates,
	 * such that repeatable results may be obtained.</p>
	 * 
	 * @param newOrderings a specification of how to shuffle the next values
	 *  to create the surrogates to generate the distribution with. The first
	 *  index is the permutation number (i.e. newOrderings.length is the number
	 *  of surrogate samples we use to bootstrap to generate the distribution here.)
	 *  Each array newOrderings[i] should be an array of length L being 
	 *  the value returned by {@link #getNumObservations()},
	 *  containing a permutation of the values in 0..(L-1).
	 * @return the distribution of AIS scores under this null hypothesis.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception where the length of each permutation in newOrderings
	 *   is not equal to the number L observations that were previously supplied.
	 */
	public EmpiricalMeasurementDistribution computeSignificance(
			int[][] newOrderings) throws Exception;

	/**
	 * Set or clear debug mode for extra debug printing to stdout
	 * 
	 * @param debug new setting for debug mode (on/off)
	 */
	public void setDebug(boolean debug);
	
	/**
	 * Return the AIS last calculated in a call to {@link #computeAverageLocalOfObservations()}
	 * or {@link #computeLocalOfPreviousObservations()} after the previous
	 * {@link #initialise()} call.
	 * 
	 * @return the last computed AIS value
	 */
	public double getLastAverage();

	/**
	 * Get the number of samples to be used for the PDFs here 
	 * which have been supplied by calls to
	 * {@link #setObservations(double[])}, {@link #addObservations(double[])}
	 * etc.
	 * 
	 * <p>Note that the number of samples is not equal to the length of time-series
	 * supplied (since we need to accumulate the first
	 * <code>(k-1)*tau + 1</code>
	 * values of each time-series before taking a sample for AIS).
	 * </p>
	 * 
	 * @return the number of samples to be used for the PDFs
	 * @throws Exception
	 */
	public int getNumObservations() throws Exception;
}
