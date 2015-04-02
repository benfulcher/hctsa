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
 * A basic interface for calculators computing measures on a univariate <i>channel</i>
 * for continuous (ie double[]) data from a
 * source to a destination time-series (ie mutual information and transfer entropy).
 * In the following, we refer to the abstract measure computed by this calculator
 * as the <i>"channel measure"</i>.
 * 
 * <p>This interface inherits from {@link ChannelCalculatorCommon} to make
 * it specific to univariate series.</p>
 * 
 * <p>Usage is as described by the paradigm for {@link ChannelCalculatorCommon},
 * but with:
 * <ul>
 * 	<li>The abstract "setObservations" and "addObservations" references therein
 *  implemented in the {@link #setObservations(double[], double[])} and
 *  {@link #addObservations(double[], double[])} etc defined here for
 *  univariate series.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @see ChannelCalculatorCommon
 */
public interface ChannelCalculator extends ChannelCalculatorCommon {

	/**
	 * Sets a single series from which to compute the PDF for the channel measure.
	 * Cannot be called in conjunction with other methods for setting/adding
	 * observations.
	 * 
	 * <p>The supplied series are certainly time-series for time-series measures
	 * such as transfer entropy, however may be simply a set of separate observations
	 * for the mutual information without a time interpretation.
	 * 
	 * @param source series of observations for the source variable. 
	 * @param destination series of observations for the destination
	 *  variable. Length must match <code>source</code>, and their indices
	 *  must correspond.
	 * @throws Exception
	 */
	public void setObservations(double source[], double destination[]) throws Exception;
	
	/**
	 * <p>Adds a new set of observations to update the PDFs with.
	 * It is intended to be called multiple times, and must
	 * be called after {@link #startAddObservations()}. Call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p><b>Important:</b> this does not append these observations to the previously
	 *  supplied observations, but treats them independently - i.e. measurements
	 *  such as the transfer entropy will not join them up to examine k
	 *  consecutive values in time.</p>
	 *  
	 * <p>Note that the arrays source and destination must not be over-written by the user
	 *  until after {@link #finaliseAddObservations()} has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param source series of observations for the source variable. 
	 * @param destination series of observations for the destination
	 *  variable. Length must match <code>source</code>, and their indices
	 *  must correspond.
	 * @throws Exception
	 */
	public void addObservations(double[] source, double[] destination) throws Exception;

	/**
	 * <p>Adds a new set of observations to update the PDFs with, as a subset
	 * of the supplied arrays.
	 * It is intended to be called multiple times, and must
	 * be called after {@link #startAddObservations()}. Call
	 * {@link #finaliseAddObservations()} once all observations have
	 * been supplied.</p>
	 * 
	 * <p><b>Important:</b> this does not append these observations to the previously
	 *  supplied observations, but treats them independently - i.e. measurements
	 *  such as the transfer entropy will not join them up to examine k
	 *  consecutive values in time.</p>
	 *  
	 * <p>Note that the arrays source and destination must not be over-written by the user
	 *  until after finaliseAddObservations() has been called
	 *  (they are not copied by this method necessarily, but the method
	 *  may simply hold a pointer to them).</p>
	 * 
	 * @param source series of observations for the source variable. 
	 * @param destination series of observations for the destination
	 *  variable. Length must match <code>source</code>, and their indices
	 *  must correspond.
	 * @param startTime first index to take observations on
	 * @param numTimeSteps number of steps from and including 
	 *   <code>startTime</code> to use
	 * @throws Exception
	 */
	public void addObservations(double[] source, double[] destination,
			int startTime, int numTimeSteps) throws Exception ;

	/**
	 * <p>Sets the single set of observations to compute the PDFs from,
	 * but only where these observations are indicated to be valid.
	 * Cannot be called in conjunction with other methods for setting/adding
	 * observations.</p>
	 * 
	 * <p>The supplied series are certainly time-series for time-series measures
	 * such as transfer entropy, however may be simply a set of separate observations
	 * for the mutual information without a time interpretation.
	 * 
	 * @param source series of observations for the source variable. 
	 * @param destination series of observations for the destination
	 *  variable. Length must match <code>source</code>, and their indices
	 *  must correspond.
	 * @param sourceValid array (with indices the same as source)
	 *  indicating whether the source at that index is valid.
	 * @param destValid array (with indices the same as destination)
	 *  indicating whether the destination at that index is valid.
	 * @throws Exception
	 */
	public void setObservations(double[] source, double[] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception;

	/**
	 * Compute the local measure values for each of the
	 * supplied samples in <code>newSourceObservations</code> and
	 * <code>newDestObservations</code>.
	 * 
	 * <p>PDFs are computed using all of the previously supplied
	 * observations, but not those in <code>newSourceObservations</code>
	 * and <code>newDestObservations</code>
	 * (unless they were some of the previously supplied samples).</p>
	 * 
	 * @param newSourceObservations series of observations for the source variable
	 *  (indexed by time or observation index)
	 * @param newDestObservations series of observations for the destination variable
	 *  (indexed by time or observation index)
	 *  Length must match <code>newSourceObservations</code>, and their indices must correspond.
	 * @return the series of local channel measure values.
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(
			double[] newSourceObservations, double[] newDestObservations)
					throws Exception;
}
