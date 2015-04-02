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
 * which may be applied to multivariate continuous time-series data.
 * That is, it is applied to <code>double[][]</code> data, indexed
 * by time then by variable number.
 * See Schreiber below for the definition of transfer entropy,
 * Lizier et al. (2011) for the extension to multivariate source
 * and destination, and
 * and Lizier et al. (2008) for the definition of local transfer entropy.
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
 * 			which at this level includes destination
 * 			history embedding length but not
 * 			delays or source embedding (TODO pull this up from TeMultivariateViaCondMi)</li>
 *		<li>Initialise the calculator using
 *			{@link #initialise()} or {@link #initialise(int, int)}
 *			or {@link #initialise(int, int, int)};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using:
 * 			<ul>
 * 				<li>{@link #setObservations(double[][], double[][])} or
 * 					{@link #setObservations(double[][], double[][], boolean[], boolean[])} or
 * 					{@link #setObservations(double[][], double[][], boolean[][], boolean[][])}
 * 					for calculations based on single time-series, OR</li>
 * 				<li>The following sequence:<ol>
 * 						<li>{@link #startAddObservations()}, then</li>
 * 						<li>One or more calls to
 * 							{@link #addObservations(double[][], double[][])} or
 * 							{@link #addObservations(double[][], double[][], int, int)}, then</li>
 * 						<li>{@link #finaliseAddObservations()};</li>
 * 					</ol></li>
 * 			</ul>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average TE: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local TE values for these samples:
 * 					{@link #computeLocalOfPreviousObservations()}</li>
 * 				<li>local TE values for a specific set of samples:
 * 				{@link #computeLocalUsingPreviousObservations(double[][], double[][])} </li>
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
 *  <li>J.T. Lizier, J. Heinzle, A. Horstmann, J.-D. Haynes, M. Prokopenko,
 *  <a href="http://dx.doi.org/10.1007/s10827-010-0271-2">
 *  "Multivariate information-theoretic measures reveal directed information
 *  structure and task relevant changes in fMRI connectivity"</a>,
 *  Journal of Computational Neuroscience, vol. 30, pp. 85-107, 2011.</li>
 * </ul>
 *
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 *
 */
public interface TransferEntropyCalculatorMultiVariate extends ChannelCalculatorMultiVariate {

	/**
	 * Property name to specify the history length k.
	 * For calculators which implement both this and
	 *  {@link TransferEntropyCalculator}, they will need
	 *  to explicitly disambiguate between {@link #K_PROP_NAME}
	 *  and {@link TransferEntropyCalculator#K_PROP_NAME}
	 */
	public static final String K_PROP_NAME = "k_HISTORY";
	
	/**
	 * Initialise the calculator for re-use with new observations.
	 * History length k, source and destination dimensions are
	 * specified here; all other parameters remain unchanged. 
	 * 
	 * @param k destination history embedding length to be considered.
	 * @param sourceDimensions number of joint variables in the source
	 * @param destDimensions number of joint variables in the destination
	 * @throws Exception
	 */
	public void initialise(int k, int sourceDimensions, int destDimensions) throws Exception;
	
}
