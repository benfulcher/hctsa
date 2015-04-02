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
 * <p>Interface for implementations of the <b>mutual information</b>,
 * which may be applied to either multivariate or merely univariate
 * continuous data.
 * That is, it is applied to <code>double[][]</code> data, where the first index
 * is observation number or time, and the second is variable number.</p>
 * 
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * 	<ol>
 * 		<li>Construct the calculator;</li>
 * 		<li>Set properties using {@link #setProperty(String, String)}
 * 			which may now include the property {@link #PROP_TIME_DIFF}
 * 			to specify a source-destination time difference (default 0);</li>
 *		<li>Initialise the calculator using {@link #initialise()} or
 *			{@link #initialise(int, int)};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using:
 * 			<ul>
 * 				<li>{@link #setObservations(double[][], double[][])},
 * 					{@link #setObservations(double[][], double[][], boolean[], boolean[])} or
 * 					{@link #setObservations(double[][], double[][], boolean[][], boolean[][])}
 * 					for calculations based on single time-series, OR</li>
 * 				<li>The following sequence:<ol>
 * 						<li>{@link #startAddObservations()}, then</li>
 * 						<li>One or more calls to {@link #addObservations(double[][], double[][])} or
 * 							{@link #addObservations(double[][], double[][], int, int)}, then</li>
 * 						<li>{@link #finaliseAddObservations()};</li>
 * 					</ol></li>
 * 			</ul>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average MI: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local MI values for these samples: {@link #computeLocalOfPreviousObservations()}</li>
 * 				<li>local MI values for a specific set of samples: {@link #computeLocalUsingPreviousObservations(double[])}</li>
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
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @see "T. M. Cover and J. A. Thomas, 'Elements of Information
Theory' (John Wiley & Sons, New York, 1991)."
 */
public interface MutualInfoCalculatorMultiVariate extends ChannelCalculatorMultiVariate {

	/**
	 * Property name for time difference between source and destination (0 by default,
	 * must be >= 0)
	 */ 
	public static final String PROP_TIME_DIFF = "TIME_DIFF";
	
	/**
	 * Compute the mutual information if the observations of the
	 * first variable (source)
	 * were ordered as per the ordering specified in newOrdering.
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
	 * @param newOrdering a specification of how to shuffle the source values
	 *  to create a surrogate source time series.
	 *  It is an array of length N (where
	 *  would be the value returned by {@link #getNumObservations()}),
	 *  containing a permutation of the values in 0..(N-1).
	 * @return the surrogate MI score if the source values were shuffled as specified.
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations(int[] newOrdering) throws Exception;
}
