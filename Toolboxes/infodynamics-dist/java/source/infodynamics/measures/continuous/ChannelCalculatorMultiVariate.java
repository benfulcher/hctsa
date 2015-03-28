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
 * A basic interface for calculators computing measures on a multivariate <i>channel</i> from a
 * source to a destination time-series (i.e. mutual information and transfer entropy).
 * In the following, we refer to the abstract measure computed by this calculator
 * as the <i>"channel measure"</i>.
 * 
 * <p>This interface inherits from {@link ChannelCalculatorCommon} to make
 * it specific to multivariate sources and destinations.</p>
 * 
 * <p>Usage is as described by the paradigm for {@link ChannelCalculatorCommon},
 * but with:
 * <ul>
 *  <li>A new {@link #initialise(int, int)} method to specify the number
 *  of joint variables in the source and destination;</li>
 * 	<li>The abstract "setObservations" and "addObservations" references therein
 *  implemented in the {@link #setObservations(double[][], double[][])} and
 *  {@link #addObservations(double[][], double[][])} etc defined here for
 *  multiivariate series.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @see ChannelCalculator
 */
public interface ChannelCalculatorMultiVariate extends ChannelCalculatorCommon {

	/**
	 * Initialise the calculator for (re-)use, with the existing
	 * (or default) values of parameters, with number of source
	 * and destination joint variables specified.
	 * Clears an PDFs of previously supplied observations.
	 *
	 * @param sourceDimensions the number of joint variables in the source
	 * @param destDimensions the number of joint variables in the destination
	 */
	public void initialise(int sourceDimensions, int destDimensions) throws Exception;

	/**
	 * Sets a single series from which to compute the PDF for the channel measure.
	 * Cannot be called in conjunction with other methods for setting/adding
	 * observations.
	 * 
	 * <p>The supplied series are certainly (multivariate) time-series for time-series measures
	 * such as transfer entropy, however may be simply a set of separate observations
	 * for the mutual information without a time interpretation.
	 * 
	 * @param source series of multivariate observations for the source variable
	 *  (first index is time or observation index, second is variable number)
	 * @param destination series of multivariate observations for the destination variable
	 *  (first index is time or observation index, second is variable number).
	 *  Length must match <code>source</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void setObservations(double[][] source, double[][] destination) throws Exception;

	/**
	 * Sets a single pair of univariate series from which to compute the PDF for the channel measure --
	 * available ONLY if both sourceDimensions and destDimensions were initialised
	 * to 1.
	 * Cannot be called in conjunction with other methods for setting/adding
	 * observations.
	 * 
	 * <p>The supplied series are certainly time-series for time-series measures
	 * such as transfer entropy, however may be simply a set of separate observations
	 * for the mutual information without a time interpretation.
	 * 
	 * <p>Design note: this method is defined here in addition to in
	 * {@link ChannelCalculator} rather than only in {@link ChannelCalculatorCommon}
	 * since {@link ConditionalTransferEntropyCalculator} etc inherit
	 * directly from {@link ChannelCalculatorCommon} but this method would not
	 * be relevant for them.
	 * </p>
	 * 
	 * @param source series of observations for the source variable. 
	 * @param destination series of observations for the destination
	 *  variable. Length must match <code>source</code>, and their indices
	 *  must correspond.
	 * @throws Exception if sourceDimensions and destDimensions were not both initialised
	 *  to 1 in {@link #initialise(int, int)}
	 */
	public void setObservations(double source[], double destination[]) throws Exception;
	
	/**
	 * <p>Adds a new set of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
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
	 * @param source series of multivariate observations for the source variable
	 *  (first index is time or observation index, second is variable number)
	 * @param destination series of multivariate observations for the destination variable
	 *  (first index is time or observation index, second is variable number).
	 *  Length must match <code>source</code>, and their indices must correspond.
	 * @throws Exception
	 */
	public void addObservations(double[][] source, double[][] destination) throws Exception;
	
	/**
	 * <p>Adds a new set of univariate observations to update the PDFs with.
	 * available ONLY if both sourceDimensions and destDimensions were initialised
	 * to 1.
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
	 * <p>Design note: this method is defined here in addition to in
	 * {@link ChannelCalculator} rather than only in {@link ChannelCalculatorCommon}
	 * since {@link ConditionalTransferEntropyCalculator} etc inherit
	 * directly from {@link ChannelCalculatorCommon} but this method would not
	 * be relevant for them.
	 * </p>
	 * 
	 * @param source series of observations for the source variable. 
	 * @param destination series of observations for the destination
	 *  variable. Length must match <code>source</code>, and their indices
	 *  must correspond.
	 * @throws Exception if sourceDimensions and destDimensions were not both initialised
	 *  to 1 in {@link #initialise(int, int)}
	 */
	public void addObservations(double[] source, double[] destination) throws Exception;

	/**
	 * <p>Adds a new sub-series of observations to update the PDFs with - is
	 * intended to be called multiple times.
	 * Must be called after {@link #startAddObservations()}; call
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
	 * @param source series of multivariate observations for the source variable
	 *  (first index is time or observation index, second is variable number)
	 * @param destination series of multivariate observations for the destination variable
	 *  (first index is time or observation index, second is variable number).
	 *  Length must match <code>source</code>, and their indices must correspond.
	 * @param startTime first index to take observations on
	 * @param numTimeSteps number of steps from and including 
	 *   <code>startTime</code> to use
	 * @throws Exception
	 */
	public void addObservations(double[][] source, double[][] destination,
			int startTime, int numTimeSteps) throws Exception;

	/**
	 * <p>Sets the single set of observations to compute the PDFs from,
	 * where the source and destination observations are valid.
	 * Cannot be called in conjunction with other methods for setting/adding
	 * observations.</p>
	 * 
	 * <p>The supplied series are certainly time-series for time-series measures
	 * such as transfer entropy, however may be simply a set of separate observations
	 * for the mutual information without a time interpretation.
	 * 
	 * @param source series of multivariate observations for the source variable
	 *  (first index is time or observation index, second is variable number)
	 * @param destination series of multivariate observations for the destination variable
	 *  (first index is time or observation index, second is variable number).
	 *  Length must match <code>source</code>, and their indices must correspond.
	 * @param sourceValid series (with first indices the same as source)
	 *  indicating whether the source at that index is valid.
	 * @param destValid time series (with first indices the same as destination)
	 *  indicating whether the destination at that index is valid.
	 */
	public void setObservations(double[][] source, double[][] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception;

	/**
	 * <p>Sets the single set of observations to compute the PDFs from,
	 * where the source and destination observations are valid for all joint variables.
	 * Cannot be called in conjunction with other methods for setting/adding
	 * observations.</p>
	 * 
	 * @param source series of multivariate observations for the source variable
	 *  (first index is time or observation index, second is variable number)
	 * @param destination series of multivariate observations for the destination variable
	 *  (first index is time or observation index, second is variable number).
	 *  Length must match <code>source</code>, and their indices must correspond.
	 * @param sourceValid series (with indices the same as source)
	 *  indicating whether each variable of the source at that index is valid.
	 * @param destValid time series (with indices the same as destination)
	 *  indicating whether each variable of the destination at that index is valid.
	 */
	public void setObservations(double[][] source, double[][] destination,
			boolean[][] sourceValid, boolean[][] destValid) throws Exception;

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
	 *  (indexed by time or observation index then variable number)
	 * @param newDestObservations series of observations for the destination variable
	 *  (indexed by time or observation index then variable number)
	 *  Length must match <code>newSourceObservations</code>, and their indices must correspond.
	 * @return the series of local channel measure values.
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(
			double[][] newSourceObservations, double[][] newDestObservations)
					throws Exception;
}
