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

package infodynamics.measures.discrete;

import infodynamics.utils.EmpiricalMeasurementDistribution;

/**
 * A basic interface for calculators computing measures on a univariate <i>channel</i>
 * for discrete (ie int[]) data from a
 * source to a destination time-series (ie mutual information and transfer entropy).
 * In the following, we refer to the abstract measure computed by this calculator
 * as the <i>"channel measure"</i>.
 * 
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * <ol>
 * 		<li>Construct the calculator;</li>
 *		<li>Initialise the calculator using {@link #initialise()} or
 *			other initialise methods defined by child classes;
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using one or more calls to
 * 			the set of {@link #addObservations(int[], int[])} methods, then</li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average channel measure: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the distribution of channel measure values under the null hypothesis
 * 					of no relationship between source and
 * 					destination values: {@link #computeSignificance(int)};</li>
 * 				<li>or other quantities as defined by child classes.</li>
 * 			</ul>
 * 		</li>
 * 		<li>
 * 		Return to step 2 to re-use the calculator on a new data set.
 * 		</li>
 * 	</ol>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public interface ChannelCalculatorDiscrete {

	/**
	 * Initialise the calculator for (re-)use, with the existing
	 * (or default) values of parameters.
	 * 
	 * @throws Exception
	 */
	public void initialise();
	
	/**
	 * <p>Adds a new set of observations to update the PDFs with.
	 * It is intended to be called multiple times.
	 * 
	 * <p><b>Important:</b> this does not append these observations to the previously
	 *  supplied observations, but treats them independently - i.e. measurements
	 *  such as the transfer entropy will not join them up to examine k
	 *  consecutive values in time.</p>
	 *  
	 * @param source series of observations for the source variable. 
	 * @param destination series of observations for the destination
	 *  variable. Length must match <code>source</code>, and their indices
	 *  must correspond.
	 */
	public void addObservations(int[] source, int[] dest);
	
	/**
	 * <p>Adds a new set of observations to update the PDFs with,
	 * from within a multivariate time-series.
	 * It is intended to be called multiple times.
	 * 
	 * <p><b>Important:</b> this does not append these observations to the previously
	 *  supplied observations, but treats them independently - i.e. measurements
	 *  such as the transfer entropy will not join them up to examine k
	 *  consecutive values in time.</p>
	 *  
	 * @param states 2D multivariate time series (first index is time
	 *  second indexes the variable)
	 * @param sourceIndex column index for the source variable. 
	 * @param destIndex column index for the destination variable.
	 */
	public void addObservations(int states[][], int sourceIndex, int destIndex);

	/**
	 * Compute the channel measure from the previously-supplied samples.
	 * 
	 * @return the estimate of the channel measure
	 */
	public double computeAverageLocalOfObservations();
	
	/**
	 * Generate a bootstrapped distribution of what the channel measure would look like,
	 * under a null hypothesis that the source values of our
	 * samples had no relation to the destination value.
	 * (Precise null hypothesis varies between MI and TE).
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for MI,
	 * conditional MI and TE.
	 * </p>
	 * 
	 * <p>Note that if several disjoint time-series have been added 
	 * as observations using {@link #addObservations(int[], int[])} etc.,
	 * then these separate "trials" will be mixed up in the generation
	 * of surrogates here.</p>
	 * 
	 * @param numPermutationsToCheck number of surrogate samples to bootstrap
	 *  to generate the distribution.
	 * @return the distribution of channel measure scores under this null hypothesis.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck);
}
