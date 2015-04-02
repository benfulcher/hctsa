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

/**
 * <p>Interface for implementations of the <b>conditional transfer entropy</b>
 * (conditional TE), which may be applied to univariate continuous
 * time-series data.
 * That is, it is applied to <code>double[]</code> data, indexed
 * by time.
 * See Schreiber below for the definition of transfer entropy,
 * and Lizier et al (2008, 2010). for the definition of local transfer entropy
 * and conditional TE, which is TE conditioned on one or 
 * more other potential sources.
 * This is also called complete TE when all other causal sources
 * are conditioned on.
 * </p>
 *  
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * 	<ol>
 * 		<li>Construct the calculator;</li>
 * 		<li>Set properties using {@link #setProperty(String, String)}
 * 			which may now include properties describing
 * 			the source and destination embedding;</li>
 *		<li>Initialise the calculator using
 *			{@link #initialise()}, {@link #initialise(int)}
 *			{@link #initialise(int, int, int)},
 *			{@link #initialise(int, int, int, int, int, int, int, int)} or
 *			{@link #initialise(int, int, int, int, int, int[], int[], int[])};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using:
 * 			<ul>
 * 				<li>{@link #setObservations(double[], double[], double[])},
 * 					{@link #setObservations(double[], double[], double[][])} or
 * 					{@link #setObservations(double[], double[], double[][], boolean[], boolean[], boolean[][])}
 * 					for calculations based on single time-series, OR</li>
 * 				<li>The following sequence:<ol>
 * 						<li>{@link #startAddObservations()}, then</li>
 * 						<li>One or more calls to
 * 							{@link #addObservations(double[], double[], double[])},
 * 							{@link #addObservations(double[], double[], double[][])},
 * 							{@link #addObservations(double[], double[], double[], int, int)} or
 * 							{@link #addObservations(double[], double[], double[][], int, int)}, then</li>
 * 						<li>{@link #finaliseAddObservations()};</li>
 * 					</ol></li>
 * 			</ul>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average conditional TE: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local TE values for these samples: {@link #computeLocalOfPreviousObservations()};</li>
 * 				<li>local TE values for a specific set of samples:
 * 				{@link #computeLocalUsingPreviousObservations(double[], double[], double[])} or
 * 				{@link #computeLocalUsingPreviousObservations(double[], double[], double[][])};</li>
 * 				<li>the distribution of MI values under the null hypothesis
 * 					of no relationship between source and
 * 					destination values: {@link #computeSignificance(int)} or
 * 					{@link #computeSignificance(int[][])}.</li>
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
 * 	<li>T. Schreiber, <a href="http://dx.doi.org/10.1103/PhysRevLett.85.461">
 * "Measuring information transfer"</a>,
 *  Physical Review Letters 85 (2) pp.461-464, 2000.</li>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href="http://dx.doi.org/10.1103/PhysRevE.77.026110">
 *  "Local information transfer as a spatiotemporal filter for complex systems"</a>
 *  Physical Review E 77, 026110, 2008.</li>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href=http://dx.doi.org/10.1063/1.3486801">
 *  "Information modification and particle collisions in distributed computation"</a>
 *  Chaos 20, 3, 037109 (2010).</li>
 * </ul>
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public interface ConditionalTransferEntropyCalculator extends ChannelCalculatorCommon {

	/**
	 * Property name to specify the history length k
	 */
	public static final String K_PROP_NAME = "k_HISTORY";	
	/**
	 * Property name for embedding delay for the destination past history vector
	 */
	public static final String K_TAU_PROP_NAME = "k_TAU";
	/**
	 * Property name for embedding length for the source past history vector
	 */
	public static final String L_PROP_NAME = "l_HISTORY";
	/**
	 * Property name for embedding delay for the source past history vector
	 */
	public static final String L_TAU_PROP_NAME = "l_TAU";
	/**
	 * Property name for source-destination delay
	 */
	public static final String DELAY_PROP_NAME = "DELAY";
	/**
	 * Property name for embedding lengths of conditional variables
	 */
	public static final String COND_EMBED_LENGTHS_PROP_NAME = "COND_EMBED_LENGTHS";
	/**
	 * Property name for embedding delays of conditional variables
	 */
	public static final String COND_EMBED_DELAYS_PROP_NAME = "COND_TAUS";
	/**
	 * Property name for conditional-destination delays of conditional variables
	 */
	public static final String COND_DELAYS_PROP_NAME = "COND_DELAYS";

	/**
	 * Initialise the calculator for re-use with new observations.
	 * A new history length k can be supplied here; all other parameters
	 * remain unchanged.
	 * 
	 * @param k destination history embedding length to be considered.
	 * @throws Exception
	 */
	public void initialise(int k) throws Exception;
	
	/**
	 * Initialise the calculator for re-use with a single conditional
	 *  variable, with the given destination,
	 *  source and conditional embedding length, setting all
	 *  embedding delays to 1, and the source-dest and
	 *  conditional-dest delays to 1.
	 *  All other parameters remain unchanged.
	 * 
	 * @param k Length of destination past history to consider
	 * @param l length of source past history to consider
	 * @param condEmbedDim embedding length for one conditional variable.
	 *  Can be 0 if there are no conditional variables.
	 * @throws Exception
	 */
	public void initialise(int k, int l, int condEmbedDim) throws Exception;

	/**
	 * Initialise the calculator with all required parameters for
	 * embeddings and delays supplied,
	 *  for a single conditional variable.
	 * All other parameters remain unchanged.
	 * 
	 * @param k Length of destination past history to consider
	 * @param k_tau embedding delay for the destination variable
	 * @param l length of source past history to consider
	 * @param l_tau embedding delay for the source variable
	 * @param delay time lag between last element of source and destination next value
	 * @param condEmbedDim embedding lengths for one conditional variable.
	 *  Can be 0 if there are no conditional variables.
	 * @param cond_tau embedding delay for the conditional variable.
	 *  Ignored if condEmbedDim == 0.
	 * @param condDelay time lags between last element of the conditional variable
	 *  and destination next value.
	 *  Ignored if condEmbedDim == 0.
	 * @throws Exception
	 */
	public void initialise(int k, int k_tau, int l, int l_tau, int delay,
			int condEmbedDim, int cond_tau, int condDelay) throws Exception;
	/**
	 * Initialise the calculator with all required parameters for
	 * embeddings and delays supplied.
	 * All other parameters remain unchanged.
	 * 
	 * @param k Length of destination past history to consider
	 * @param k_tau embedding delay for the destination variable
	 * @param l length of source past history to consider
	 * @param l_tau embedding delay for the source variable
	 * @param delay time lag between last element of source and destination next value
	 * @param condEmbedDims array of embedding lengths for each conditional variable.
	 *  Can be an empty array or null if there are no conditional variables.
	 * @param cond_taus array of embedding delays for the conditional variables.
	 *  Must be same length as condEmbedDims array.
	 * @param condDelays array of time lags between last element of each conditional variable
	 *  and destination next value.
	 *  Must be same length as condEmbedDims array.
	 * @throws Exception for inconsistent arguments, e.g. if array lengths differ between 
	 *  condEmbedDims, cond_taus and condDelays.
	 */
	public void initialise(int k, int k_tau, int l, int l_tau, int delay,
			int[] condEmbedDims, int[] cond_taus, int[] condDelays) throws Exception;
	
	// Overriding the javadocs here, the method is already defined on ChannelCalculatorCommon
	/**
	 * Sets properties for the conditional TE calculator.
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 *  
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 * 		<li>{@link #K_PROP_NAME} -- destination history embedding length k
	 * 			(default value 1)</li>
	 * 		<li>{@link #K_TAU_PROP_NAME} -- embedding delay for the destination past history vector
	 * 			(default value 1)</li>
	 * 		<li>{@link #L_PROP_NAME} -- embedding length for the source past history vector
	 * 			(default value 1)</li>
	 * 		<li>{@link #L_TAU_PROP_NAME} -- embedding delay for the source past history vector
	 * 			(default value 1)</li>
	 * 		<li>{@link #DELAY_PROP_NAME} -- source-destination delay (default value is 1)</li>
	 * 		<li>{@link #COND_EMBED_LENGTHS_PROP_NAME} -- conditional variables embedding lengths
	 * 			as a comma separated integer list (default value "1" for a single variable)</li>
	 * 		<li>{@link #COND_EMBED_DELAYS_PROP_NAME} -- conditional variables embedding delays
	 * 			as a comma separated integer list (default value "1" for a single variable)</li>
	 * 		<li>{@link #COND_DELAYS_PROP_NAME} --  delays from the conditional variables to the destination
	 * 			as a comma separated integer list (default value "1" for a single variable)</li>
	 * </ul>
	 * 
	 * <p>While {@link #COND_EMBED_LENGTHS_PROP_NAME},
	 * {@link #COND_EMBED_DELAYS_PROP_NAME} and {@link #COND_DELAYS_PROP_NAME}
	 * may be set separately and may therefore be arrays of different
	 * lengths, {@link #COND_EMBED_LENGTHS_PROP_NAME} acts as the master
	 * from which we determine the number of conditional variables that we have,
	 * and we will check that they have the same lengths at the next initialisation.</p>
	 * 
	 * <p><b>Note:</b> further properties may be defined by child classes.</p>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property.
	 * @throws Exception if there is a problem with the supplied value, 
	 * or if the property is recognised but unsupported.
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue) throws Exception;
	
	/**
	 * Sets the single time-series from which to compute the PDFs.
	 * Cannot be called in conjunction with other methods for setting/adding
	 * observations.
	 * 
	 * @param source time-series of observations for the source variable. 
	 * @param destination time-series of observations for the destination
	 *  variable. Length must match <code>source</code>.
	 * @param conditionals 2D time-series array for the conditional variables
	 *        (first index is time, second index is variable number).
	 *        Length must match <code>source</code>.
	 * @throws Exception
	 */
	public void setObservations(double[] source, double[] destination, double[][] conditionals) throws Exception;
	
	/**
	 * Sets the single time-series from which to compute the PDFs.
	 * Cannot be called in conjunction with other methods for setting/adding
	 * observations.
	 * 
	 * @param source time-series of observations for the source variable. 
	 * @param destination time-series of observations for the destination
	 *  variable. Length must match <code>source</code>.
	 * @param conditionals 1D time series array for the conditional variables
	 *        (indexed by time only) -- valid only if the calculator
	 *        was initialised for a single conditional variable.
	 *        Length must match <code>source</code>.
	 * @throws Exception for example if the calculator was not initialised for
	 *        a single conditional variable
	 */
	public void setObservations(double[] source, double[] destination, double[] conditionals) throws Exception;

	/**
	 * <p>Adds a new set of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p><b>Important:</b> this does not append/concatenate these observations to the previously
	 *  supplied observations, but treats them independently.</p>
	 *  
	 * <p>Note that the arrays source, destination and conditionals must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param source time-series observations for the source variable
	 * @param destination time-series of observations for the destination
	 *  variable. Length must match <code>source</code>.
	 * @param conditionals 2D time-series array for the conditional variables
	 *        (first index is time, second index is variable number).
	 *        Length must match <code>source</code>.
	 * @throws Exception
	 */
	public void addObservations(double[] source, double[] destination,  double[][] conditionals) throws Exception;

	/**
	 * <p>Adds a new set of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p><b>Important:</b> this does not append/concatenate these observations to the previously
	 *  supplied observations, but treats them independently.</p>
	 *  
	 * <p>Note that the arrays source, destination and conditionals must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param source time-series observations for the source variable
	 * @param destination time-series of observations for the destination
	 *  variable. Length must match <code>source</code>.
	 * @param conditionals 1D time series array for the conditional variables
	 *        (indexed by time only) -- valid only if the calculator
	 *        was initialised for a single conditional variable.
	 *        Length must match <code>source</code>.
	 * @throws Exception for example if the calculator was not initialised for
	 *        a single conditional variable
	 */
	public void addObservations(double[] source, double[] destination,  double[] conditionals) throws Exception;

	/**
	 * <p>Adds a new sub-series of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p><b>Important:</b> this does not append/concatenate these observations to the previously
	 *  supplied observations, but treats them independently.</p>
	 *  
	 * <p>Note that the arrays source, destination and conditionals must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param source time-series observations for the source variable
	 * @param destination time-series of observations for the destination
	 *  variable. Length must match <code>source</code>.
	 * @param conditionals 2D time-series array for the conditional variables
	 *        (first index is time, second index is variable number).
	 *        Length must match <code>source</code>.
	 * @param startTime first time index to take observations on
	 * @param numTimeSteps number of time steps to use
	 * @throws Exception
	 */
	public void addObservations(double[] source, double[] destination,
			double[][] conditionals,
			int startTime, int numTimeSteps) throws Exception ;

	/**
	 * <p>Adds a new sub-series of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p><b>Important:</b> this does not append/concatenate these observations to the previously
	 *  supplied observations, but treats them independently.</p>
	 *  
	 * <p>Note that the arrays source, destination and conditionals must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param source time-series observations for the source variable
	 * @param destination time-series of observations for the destination
	 *  variable. Length must match <code>source</code>.
	 * @param conditionals 1D time-series array for the conditional variables
	 *        (indexed by time only) -- valid only if the calculator
	 *        was initialised for a single conditional variable.
	 *        Length must match <code>source</code>.
	 * @param startTime first time index to take observations on
	 * @param numTimeSteps number of time steps to use
	 * @throws Exception for example if the calculator was not initialised for
	 *        a single conditional variable
	 */
	public void addObservations(double[] source, double[] destination,
			double[] conditionals,
			int startTime, int numTimeSteps) throws Exception ;

	/**
	 * <p>Sets the single set of observations to compute the PDFs from,
	 * but only where these observations are indicated to be valid.
	 * Cannot be called in conjunction with 
	 * {@link #startAddObservations()}/{@link #addObservations(double[], double[])} /
	 * {@link #finaliseAddObservations()}.</p>
	 * 
	 * @param source time-series observations for the source variable
	 * @param destination time-series of observations for the destination
	 *  variable. Length must match <code>source</code>.
	 * @param conditionals 2D time-series array for the conditional variables
	 *        (first index is time, second index is variable number).
	 *        Length must match <code>source</code>.
	 * @param sourceValid time-series (with time indices the same as source)
	 *  indicating whether the source at that point is valid.
	 * @param destValid time-series (with time indices the same as destination)
	 *  indicating whether the destination at that point is valid.
	 * @param conditionalsValid 2D time-series (with time indices the same as conditionals)
	 *  indicating whether the conditional variables at that point are valid.
	 * @throws Exception
	 */
	public void setObservations(double[] source, double[] destination,
			double[][] conditionals,
			boolean[] sourceValid, boolean[] destValid, boolean[][] conditionalsValid) throws Exception;

	/**
	 * Compute local conditional transfer entropy values for the
	 *  observations in the given time-series,
	 *  using the PDFs computed from the previously supplied method calls.
	 * 
	 * @param newSourceObservations new time-series observations for the source variable
	 * @param newDestObservations new time-series of observations for the destination
	 *  variable. Length must match <code>source</code>.
	 * @param newCondObservations new 2D time-series array for the conditional variables
	 *        (first index is time, second index is variable number).
	 *        Length must match <code>source</code>.
	 * @return the time-series of local conditional TE values.
	 *        First values will be set to zero until enough samples are
	 *        accumulated to embed the variables as per the parameter settings.
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(
			double[] newSourceObservations, double[] newDestObservations,
			double[][] newCondObservations) throws Exception;

	/**
	 * Compute local conditional transfer entropy values for the
	 *  observations in the given parameters,
	 *  using the PDFs computed from the previously supplied method calls.
	 * 
	 * @param newSourceObservations new time-series observations for the source variable
	 * @param newDestObservations new time-series of observations for the destination
	 *  variable. Length must match <code>source</code>.
	 * @param newCondObservations new 1D time-series array for the conditional variables
	 *        (indexed by time only) -- valid only if the calculator
	 *        was initialised for a single conditional variable.
	 *        Length must match <code>source</code>.
	 * @return the time-series of local conditional TE values.
	 *        First values will be set to zero until enough samples are
	 *        accumulated to embed the variables as per the parameter settings.
	 * @throws Exception for example if the calculator was not initialised for
	 *        a single conditional variable
	 */
	public double[] computeLocalUsingPreviousObservations(
			double[] newSourceObservations, double[] newDestObservations,
			double[] newCondObservations) throws Exception;
}
