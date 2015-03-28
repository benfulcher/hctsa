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
 * Interface for implementations of
 * entropy estimators on continuous multivariate data
 * (ie <code>double[][]</code> arrays
 * where first index is time or observation number, second is variable index).
 * 
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * 	<ol>
 * 		<li>Construct the calculator;</li>
 * 		<li>Set properties using {@link #setProperty(String, String)};</li>
 *		<li>Initialise the calculator {@link #initialise(int)};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      to set up the PDFs, using:
 * 			{@link #setObservations(double[][])};</li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average entropy: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local entropy values for these samples: {@link #computeLocalOfPreviousObservations()}</li>
 * 				<li>local entropy values for a specific set of samples:
 * 					{@link #computeLocalUsingPreviousObservations(double[])}.</li>
 * 			</ul>
 * 		</li>
 * 		<li>Return to step 2 or 3 to re-use the calculator on a new data set.</li>
 * 	</ol>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. M. Cover and J. A. Thomas, 'Elements of Information
Theory' (John Wiley & Sons, New York, 1991).</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public interface EntropyCalculatorMultiVariate {

	/**
	 * Initialise the calculator for (re-)use, with the existing (or default) values
	 * of calculator-specific parameters.
	 * Clears any PDFs of previously supplied observations.
	 * 
	 * @param dimensions number of joint variables to be investigated
	 */
	public void initialise(int dimensions);
	
	/**
	 * Set properties for the underlying calculator implementation.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * <p>No general properties are defined at the interface level, i.e.
	 * there are only calculator-specific properties.</p>
	 *  
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception;

	/**
	 * Set the observations for which to compute the PDFs for the entropy
	 * Should only be called once, the last call contains the
	 *  observations that are used (they are not accumulated).
	 * 
	 * @param observations multivariate time series of observations; first index
	 *  is time step, second index is variable number (total should match dimensions
	 *  supplied to {@link #initialise(int)}
	 * @throws Exception if the dimensions of the observations do not match 
	 *  the expected value supplied in {@link #initialise(int)}; implementations
	 *  may throw other more specific exceptions also.
	 */
	public void setObservations(double observations[][]) throws Exception;
	
	/**
	 * Compute the entropy from the previously-supplied samples.
	 * 
	 * @return the entropy estimate, in bits or nats depending on the estimator.
	 */
	public double computeAverageLocalOfObservations();
	
	/**
	 * Compute the local entropy values for each of the
	 * supplied samples in <code>newObservations</code>.
	 * 
	 * <p>PDFs are computed using all of the previously supplied
	 * observations, but not those in <code>newObservations</code> (unless they were
	 * some of the previously supplied samples).</p>
	 * 
	 * @param newObservations multivariate time-series for which to compute
	 * 	local entropy values (see {@link #setObservations(double[][])} for its format)
	 * @return time-series of local entropy values corresponding to each entry
	 * in <code>newObservations</code>
	 */
	public double[] computeLocalUsingPreviousObservations(double newObservations[][]) throws Exception;

	/**
	 * Compute the local entropy values for each of the
	 * previously-supplied samples.
	 * 
	 * <p>PDFs are computed using all of the previously supplied
	 * observations.</p>
	 * 
	 * @return the array of local entropy values corresponding to each
	 *  previously supplied observation.
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception;
	
	/**
	 * Return the entropy last calculated in a call to {@link #computeAverageLocalOfObservations()}
	 * or {@link #computeLocalOfPreviousObservations()} after the previous
	 * {@link #initialise(int)} call.
	 * 
	 * @return the last computed entropy value
	 */
	public double getLastAverage();
	
	/**
	 * Set or clear debug mode for extra debug printing to stdout
	 * 
	 * @param debug new setting for debug mode (on/off)
	 */
	public void setDebug(boolean debug);
	
	/**
	 * Get the number of samples to be used for the PDFs here 
	 * which have been supplied by calls to
	 * {@link #setObservations(double[][])}.
	 * 
	 * @return the number of samples to be used for the PDFs
	 * @throws Exception
	 */
	public int getNumObservations() throws Exception;
}
