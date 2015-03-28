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
 * <p>Interface for implementations of the multi-information or integration,
 * which may be applied to either multivariate or merely univariate
 * continuous data.
 * That is, it is applied to <code>double[][]</code> data, where the first index
 * is observation number or time, and the second is variable number.
 * See Tononi et al. below for the definition of multi-information/integration.</p>
 * 
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * 	<ol>
 * 		<li>Construct the calculator;</li>
 * 		<li>Set properties using {@link #setProperty(String, String)};</li>
 *		<li>Initialise the calculator using {@link #initialise(int)};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using:
 * 			<ul>
 * 				<li>{@link #setObservations(double[][])}
 * 					for calculations based on single time-series, OR</li>
 * 				<li>The following sequence:<ol>
 * 						<li>{@link #startAddObservations()}, then</li>
 * 						<li>One or more calls to {@link #addObservations(double[][])} or
 * 							{@link #addObservation(double[])}, then</li>
 * 						<li>{@link #finaliseAddObservations()};</li>
 * 					</ol></li>
 * 			</ul>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average multi-information:
 *					{@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local multi-information values for these samples:
 * 					{@link #computeLocalOfPreviousObservations()}</li>
 * 				<li>local multi-information values for a specific set of samples:
 * 					{@link #computeLocalUsingPreviousObservations(double[][])}</li>
 * 			</ul>
 * 		</li>
 * 		<li>
 * 		Return to step 2 or 3 to re-use the calculator on a new data set.
 * 		</li>
 * 	</ol>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>"G. Tononi, O. Sporns, G. M. Edelman,
 *  <a href="http://dx.doi.org/10.1073/pnas.91.11.5033">"A measure for
 *  brain complexity:
 * 	relating functional segregation and integration in the nervous system"</a>
 *  Proceedings of the National Academy of Sciences, Vol. 91, No. 11.
 *  (1994), pp. 5033-5037.</li>
 * </ul>
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public interface MultiInfoCalculator {

	/**
	 * Property name for whether to normalise incoming values to mean 0,
	 * standard deviation 1 (default true)
	 */
	public static final String PROP_NORMALISE = "NORMALISE";
	/**
	 * Property name for whether to use less than 100% of the samples
	 * in making the calculation; double value for this property in 0..1
	 * gives the proportion of values to use (default is 1 meaning all).
	 */
	public static final String SAMPLING_FACTOR_PROP_NAME = "SAMPLING_FACTOR";

	/**
	 * Initialise the calculator for (re-)use, with the existing
	 * (or default) values of parameters, with number of 
	 * joint variables specified.
	 * Clears an PDFs of previously supplied observations.
	 *
	 * @param dimensions the number of joint variables to consider
	 */
	public void initialise(int dimensions);

	/**
	 * Set properties for the calculator.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 * 		<li>{@link #PROP_NORMALISE} -- whether to normalise the incoming variables 
	 * 			to mean 0, standard deviation 1, or not (default false).</li>
	 * 		<li>{@link #SAMPLING_FACTOR_PROP_NAME} -- whether to use less than 100% of the samples
	 * 			in making the calculation; double value for this property in 0..1
	 * 			gives the proportion of values to use (default is 1 meaning all).</li>
	 * 		<li>Any properties defined by the relevant child interfaces and classes.</li>
	 * </ul>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception;
	
	/**
	 * Sets a single series from which to compute the PDF for the multi-information.
	 * Cannot be called in conjunction with other methods for setting/adding
	 * observations.
	 * 
	 * <p>The supplied series may be a time-series, or may be simply
	 * a set of separate observations
	 * without a time interpretation.</p>
	 * 
	 * <p>Should only be called once, the last call contains the
	 *  observations that are used (they are not accumulated).</p> 
	 * 
	 * @param observations series of multivariate observations
	 *  (first index is time or observation index, second is variable number)
	 * @throws Exception
	 */
	public void setObservations(double observations[][]) throws Exception;
	
	/**
	 * Signal that we will add in the samples for computing the PDF 
	 * from several disjoint time-series or trials via calls to
	 * {@link #addObservation(double[])} or {@link #addObservations(double[][])}
	 * rather than {@link #setDebug(boolean)}.
	 */
	public void startAddObservations();
	
	/**
	 * <p>Adds a new (single) observation to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p>Note that the arrays must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param observation a single multivariate observation
	 *  (index is variable number)
	 */
	public void addObservation(double observation[]);
	
	/**
	 * <p>Adds a new set of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p>Note that the arrays must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param observations series of multivariate observations
	 *  (first index is time or observation index, second is variable number)
	 */
	public void addObservations(double observations[][]);
	
	/**
	 * Signal that the observations are now all added, PDFs can now be constructed.
	 * 
	 * @throws Exception 
	 */
	public void finaliseAddObservations() throws Exception;

	/**
	 * Compute the multi-information from the previously-supplied samples.
	 * 
	 * @return the estimate of the multi-information
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations() throws Exception;

	/**
	 * <p>Computes the local values of the multi-information,
	 *  for each valid observation in the previously supplied observations
	 *  (with PDFs computed using all of the previously supplied observation sets).</p>
	 *  
	 * <p>If the samples were supplied via a single call i.e.
	 * {@link #setObservations(double[][])},
	 * then the return value is a single <i>time-series</i> of local
	 * channel measure values corresponding to these samples.</p>
	 * 
	 * <p>Otherwise where disjoint time-series observations were supplied using several 
	 *  calls such as {@link addObservations(double[][])}
	 *  then the local values for each disjoint observation set will be appended here
	 *  to create a single "time-series" return array.</p>
	 *  
	 * @return the "time-series" of local multi-information values.
	 * @throws Exception
	 */
	public double[] computeLocalOfPreviousObservations()
		throws Exception;
	
	/**
	 * Compute the local multi-information values for each of the
	 * supplied samples in <code>states</code>.
	 * 
	 * <p>PDFs are computed using all of the previously supplied
	 * observations, but not those in <code>states</code>
	 * (unless they were
	 * some of the previously supplied samples).</p>
	 * 
	 * @param states series of multivariate observations
	 *  (first index is time or observation index, second is variable number)
	 * @return the series of local multi-information values.
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(double states[][])
		throws Exception;

	/**
	 * Compute the multi information would look like were all time series
	 *  (bar the first) reordered
	 *  as per the array of time indices in newOrdering.
	 * 
	 * <p>The reordering array contains the reordering for each marginal variable (first index).
	 * The user should ensure that all values 0..N-1 are represented exactly once in the
	 *  array reordering and that no other values are included here.</p>
	 *  
	 * <p>Note that if several disjoint time-series have been added 
	 * as observations using {@link #addObservations(double[])} etc.,
	 * then these separate "trials" will be mixed up in the generation
	 * of a shuffled source series here.</p>
	 * 
	 * <p>This method is primarily intended for use in {@link #computeSignificance(int[][])}
	 * however has been made public in case users wish to access it.
	 * </p>
	 * 
	 * @param newOrdering the specific permuted new orderings to use. First index is the variable number
	 *  (minus 1, since we don't reorder the first variable),
	 *  second index is the time step, the value is the reordered time step to use
	 *  for that variable at the given time step.
	 *  The values must be an array of length N (where
	 *  would be the value returned by {@link #getNumObservations()}),
	 *  containing a permutation of the values in 0..(N-1).
	 *  If null, no reordering is performed.
	 * @return what the average multi-info would look like under this reordering
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations(int[][] newOrdering) throws Exception;

	/**
	 * Generate a resampled distribution of what the multi-information would look like,
	 * under a null hypothesis that the individual values of each
	 * variable in the 
	 * samples have no relation to eachother.
	 * That is, we destroy the p(x,y,z,..) correlations, while
	 * retaining the p(x), p(y),.. marginals, to check how
	 *  significant this multi-information actually was.
	 *  
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for MI,
	 * we are extending that here.
	 * </p>
	 * 
	 * <p>Note that if several disjoint time-series have been added 
	 * as observations using {@link #addObservations(double[])} etc.,
	 * then these separate "trials" will be mixed up in the generation
	 * of surrogates here.</p>
	 * 
	 * <p>This method (in contrast to {@link #computeSignificance(int[][][])})
	 * creates <i>random</i> shufflings of the next values for the surrogate MultiInfo
	 * calculations.</p>
	 * 
	 * @param numPermutationsToCheck number of surrogate samples to permute
	 *  to generate the distribution.
	 * @return the distribution of surrogate multi-info values under this null hypothesis.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception;
	
	/**
	 * Generate a resampled distribution of what the multi-information would look like,
	 * under a null hypothesis that the individual values of each
	 * variable in the 
	 * samples have no relation to eachother.
	 * That is, we destroy the p(x,y,z,..) correlations, while
	 * retaining the p(x), p(y),.. marginals, to check how
	 *  significant this multi-information actually was.
	 *  
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for MI,
	 * we are extending that here.
	 * </p>
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
	 * @param newOrderings a specification of how to shuffle the values
	 *  to create the surrogates to generate the distribution with. The first
	 *  index is the permutation number (i.e. newOrderings.length is the number
	 *  of surrogate samples we use to bootstrap to generate the distribution here.)
	 *  The second index is the variable number (minus 1, since we don't reorder
	 *  the first variable),
	 *  Each array newOrderings[i][v] should be an array of length N (where
	 *  would be the value returned by {@link #getNumObservations()}),
	 *  containing a permutation of the values in 0..(N-1).
	 * @return the distribution of surrogate multi-info values under this null hypothesis.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception where the length of each permutation in newOrderings
	 *   is not equal to the number N samples that were previously supplied.
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][][] newOrderings) throws Exception;
	
	/**
	 * Set or clear debug mode for extra debug printing to stdout
	 * 
	 * @param debug new setting for debug mode (on/off)
	 */
	public void setDebug(boolean debug);
	
	/**
	 * Return the multi-information last calculated in a call to {@link #computeAverageLocalOfObservations()}
	 * or {@link #computeLocalOfPreviousObservations()} after the previous
	 * {@link #initialise(int)} call.
	 * 
	 * @return the last computed multi-information value
	 */
	public double getLastAverage();
}
