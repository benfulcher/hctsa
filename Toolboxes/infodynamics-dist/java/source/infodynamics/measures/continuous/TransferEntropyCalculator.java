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
 * <p>Interface for implementations of the <b>transfer entropy</b> (TE),
 * which may be applied to univariate continuous time-series data.
 * That is, it is applied to <code>double[]</code> data, indexed
 * by time.
 * See Schreiber below for the definition of transfer entropy,
 * and Lizier et al. for the definition of local transfer entropy.
 * Specifically, this class implements the pairwise or <i>apparent</i>
 * transfer entropy; i.e. we compute the transfer that appears to
 *  come from a single source variable, without examining any other
 *  potential sources
 *  (see Lizier et al, PRE, 2008).</p>
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
 *			{@link #initialise()} or {@link #initialise(int)};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using:
 * 			<ul>
 * 				<li>{@link #setObservations(double[], double[])} or
 * 					{@link #setObservations(double[], double[], boolean[], boolean[])}
 * 					for calculations based on single time-series, OR</li>
 * 				<li>The following sequence:<ol>
 * 						<li>{@link #startAddObservations()}, then</li>
 * 						<li>One or more calls to
 * 							{@link #addObservations(double[], double[])} or
 * 							{@link #addObservations(double[], double[], int, int)}, then</li>
 * 						<li>{@link #finaliseAddObservations()};</li>
 * 					</ol></li>
 * 			</ul></li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average TE: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local TE values for these samples: {@link #computeLocalOfPreviousObservations()}</li>
 * 				<li>local TE values for a specific set of samples:
 * 				{@link #computeLocalUsingPreviousObservations(double[], double[])} </li>
 * 				<li>the distribution of TE values under the null hypothesis
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
 * </ul>
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public interface TransferEntropyCalculator extends ChannelCalculator {

	/**
	 * Property name to specify the destination history embedding length k
	 * (default value 1)
	 */
	public static final String K_PROP_NAME = "k_HISTORY";
	/**
	 * Property name for embedding delay for the destination past history vector
	 * (default value 1)
	 */
	public static final String K_TAU_PROP_NAME = "k_TAU";
	/**
	 * Property name for embedding length for the source past history vector
	 * (default value 1)
	 */
	public static final String L_PROP_NAME = "l_HISTORY";
	/**
	 * Property name for embedding delay for the source past history vector
	 * (default value 1)
	 */
	public static final String L_TAU_PROP_NAME = "l_TAU";
	/**
	 * Property name for source-destination delay (default value is 1)
	 */
	public static final String DELAY_PROP_NAME = "DELAY";
	
	/**
	 * Initialise the calculator for re-use with new observations.
	 * A new history length k can be supplied here; all other parameters
	 * remain unchanged.
	 * 
	 * @param k destination history embedding length to be considered.
	 * @throws Exception
	 */
	public void initialise(int k) throws Exception;
	
	// Overriding the javadocs here, the method is already defined on ChannelCalculatorCommon
	/**
	 * Sets properties for the TE calculator.
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
	 * </ul>
	 * <p><b>Note:</b> further properties may be defined by child classes.</p>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property.
	 * @throws Exception if there is a problem with the supplied value, 
	 * or if the property is recognised but unsupported (eg some
	 * calculators do not support all of the embedding properties).
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue) throws Exception;

	// Overriding the javadocs here, the method is already defined on ChannelCalculatorCommon
	/**
	 * @return the "time-series" of local TE values.
	 *  If the samples were supplied via a single call such as
	 *  {@link #setObservations(double[])},
	 *  then this time-series will have zero entries until the source
	 *  and destination embedding vectors are defined.
	 */
	@Override
	public double[] computeLocalOfPreviousObservations() throws Exception;
}
