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

/**
 * <p>Base class for our information-theoretic calculators
 *  on discrete (int[]) data,
 *  providing common functionality
 *  for user-level measure classes.</p>
 * 
 * <p>
 * Usage of the child classes extending this class is intended to follow this paradigm:
 * </p>
 * <ol>
 * 		<li>Construct the calculator;</li>
 *		<li>Initialise the calculator using {@link #initialise()} or
 *			other initialise methods defined by child classes;</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using one or more calls to
 * 			sets of "addObservations" methods defined by child classes, then</li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average measure: {@link #computeAverageLocalOfObservations()};</li>
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
public abstract class InfoMeasureCalculatorDiscrete {

	/**
	 * Last computed average of the measure
	 */
	protected double average = 0.0;
	/**
	 * Last computed max local value of the measure
	 */
	protected double max = 0.0;
	/**
	 * Last computed min local value of the measure
	 */
	protected double min = 0.0;
	/**
	 * Last computed standard deviation of local values of the measure
	 */
	protected double std = 0.0;
	/**
	 * Number of observations supplied for the PDFs
	 */
	protected int observations = 0;
	/**
	 * Number of available quantised states for each variable
	 * (ie binary is base-2).
	 */
	protected int base = 0; // number of individual states. Need initialised to 0 for changedSizes

	/**
	 * Cached value of ln(base)
	 */
	protected double log_base = 0;
	/**
	 * Cached value of ln(2)
	 */
	protected double log_2 = Math.log(2.0);
	/**
	 * Cache of whether the base is a power of 2
	 */
	protected boolean power_of_2_base = false;
	/**
	 * Cached value of log_2(base)
	 */
	protected int log_2_base = 0;
	/**
	 * Whether we're in debug mode
	 */
	protected boolean debug = false;
	
	/**
	 * Construct an instance
	 * 
	 * @param base number of quantisation levels for each variable.
	 *        E.g. binary variables are in base-2.
	 */
	protected InfoMeasureCalculatorDiscrete(int base) {
		
		this.base = base;
		log_base = Math.log(base);

		if (base < 2) {
			throw new RuntimeException("Can't calculate info theoretic measures for base " + base);
		}
		
		// Check if we've got a power of 2
		power_of_2_base = isPowerOf2(base);
		if (power_of_2_base) {
			log_2_base = (int) Math.round(Math.log(base) / Math.log(2));
		}
	}
	
	/**
	 * Initialise the calculator for re-use with new observations.
	 * (Child classes should clear the existing PDFs)
	 */
	public void initialise(){
		average = 0.0;
		max = 0.0;
		min = 0.0;
		std = 0.0;
		observations = 0;
	}
	
	/**
	 * Return the measure last calculated in a call to
	 * {@link #computeAverageLocalOfObservations()}
	 * or related methods after the previous
	 * {@link #initialise()} call.
	 * 
	 * @return the last computed measure value
	 */
	public final double getLastAverage() {
		return average;
	}

	/**
	 * Return the last computed max local value of the measure.
	 * Not declaring this final so that separable calculator
	 *  can throw an exception on it since it does not support it
	 */
	public double getLastMax() {
		return max;
	}

	/**
	 * Return the last computed min local value of the measure.
	 * Not declaring this final so that separable calculator
	 *  can throw an exception on it since it does not support it
	 */
	public double getLastMin() {
		return min;
	}
	
	/**
	 * Return the last computed standard deviation of
	 * local values of the measure.
	 */
	public final double getLastStd() {
		return std;
	}
	
	/**
	 * Get the number of samples to be used for the PDFs here 
	 * which have been supplied by calls to
	 * "setObservations", "addObservations" etc.
	 * 	
	 * <p>Note that the number of samples may not be equal to the length of time-series
	 * supplied (e.g. for transfer entropy, where we need to accumulate
	 * a number of samples for the past history of the destination).
	 * </p>
	 * 
	 * @return the number of samples to be used for the PDFs
	 */
	public final int getNumObservations() {
		return observations;
	}
	
	/**
	 * Return whether a given integer is a power of 2
	 * 
	 * @param num an integer
	 * @return whether the integer is a power of 2
	 */
	public final static boolean isPowerOf2(int num) {
		int bits = 0;
		int shiftedValue = num;
		for (int b = 0; b < Integer.SIZE; b ++) {
			if ((shiftedValue & 0x01) > 0) {
				// LSB is a 1
				bits++;
				if (bits > 1) {
					// Encountered more than 1 bit set to 1.
					// Is not a power of 2
					return false;
				}
			}
			// Shift a new bit down into the LSB
			shiftedValue = shiftedValue >> 1;
		}
		// Post: num has either 1 bit set to 1 or 0 bits set to 1
		if (bits == 0) {
			return false;
		}
		return true;
	}
	
	/**
	 * Compute the average value of the measure
	 * from the previously-supplied samples.
	 * 
	 * @return the estimate of the measure
	 */
	public abstract double computeAverageLocalOfObservations();
	
	/**
	 * Set or clear debug mode for extra debug printing to stdout
	 * 
	 * @param debug new setting for debug mode (on/off)
	 */
	public void setDebug(boolean debug) {
		this.debug = debug;
	}
}
