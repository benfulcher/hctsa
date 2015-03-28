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
 * entropy estimators on continuous univariate data
 * (ie <code>double[]</code> arrays).
 *  
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * 	<ol>
 * 		<li>Construct the calculator;</li>
 * 		<li>Set properties using {@link #setProperty(String, String)};</li>
 *		<li>Initialise the calculator {@link #initialise()};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      to set up the PDFs, using: {@link #setObservations(double[])};</li>
 * 		<li>Compute the entropy using:
 * 			{@link #computeAverageLocalOfObservations()};</li>
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
public interface EntropyCalculator {

	/**
	 * Initialise the calculator for (re-)use, with the existing (or default) values of calculator-specific parameters
 	 * Clears any PDFs of previously supplied observations.
	 */
	public void initialise() throws Exception;
	
	/**
	 * Set properties for the underlying calculator implementation.
	 * New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * <p>No general properties are defined at the interface level, i.e.
	 * there are only calculator-specific properties.</p>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	public void setProperty(String propertyName, String propertyValue) throws Exception;

	// TODO Add addObservations() methods for entropy calculator
	
	/**
	 * Sets the samples from which to compute the PDF for the entropy.
	 * Should only be called once, the last call contains the
	 *  observations that are used (they are not accumulated). 
	 * 
	 * @param observations array of (univariate) samples
	 */
	public void setObservations(double observations[]);
	
	/**
	 * Compute the entropy from the previously-supplied samples.
	 * 
	 * @return the entropy estimate
	 */
	public double computeAverageLocalOfObservations();

	/**
	 * Set or clear debug mode for extra debug printing to stdout
	 * 
	 * @param debug new setting for debug mode (on/off)
	 */
	public void setDebug(boolean debug);
}
