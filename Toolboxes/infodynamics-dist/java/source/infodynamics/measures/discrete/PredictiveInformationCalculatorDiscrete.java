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

import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

/**
 * <p> Predictive information calculator for univariate discrete (int[]) data.
 * See definition of Predictive information (PI) by Bialek et al. below,
 * also is a form of the Excess Entropy (see Crutchfield and Feldman below),
 * Basically, PI is the mutual information between the past <i>state</i>
 * of a time-series process <i>X</i> and its future <i>state</i>.
 * The past <i>state</i> at time <code>n</code>
 * is represented by an embedding vector of <code>k</code> values from <code>X_n</code> backwards,
 * each separated by <code>\tau</code> steps, giving
 * <code><b>X^k_n</b> = [ X_{n-(k-1)\tau}, ... , X_{n-\tau}, X_n]</code>.
 * We call <code>k</code> the embedding dimension, and <code>\tau</code>
 * the embedding delay (only delay = 1 is implemented at the moment).
 * The future <i>state</i> at time <code>n</code>
 * is defined similarly into the future:
 * each separated by <code>\tau</code> steps, giving
 * <code><b>X^k+_n</b> = [ X_{n+1}, X_{n+\tau}, X_{n+(k-1)\tau}]</code>.
 * PI is then the mutual information between <b>X^k_n</b> and <b>X^k+_n</b>.</p>
 * 
 * <p>Usage of the class is intended to follow this paradigm:</p>
 * <ol>
 * 		<li>Construct the calculator:
 * 			{@link #PredictiveInformationCalculatorDiscrete(int, int)};</li>
 *		<li>Initialise the calculator using {@link #initialise()};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using one or more calls to
 * 			sets of {@link #addObservations(int[])} methods, then</li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average entropy: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>local entropy values, such as {@link #computeLocal(int[])};</li>
 * 				<li>and variants of these.</li>
 * 			</ul>
 * 		</li>
 * 		<li>As an alternative to steps 3 and 4, the user may undertake
 * 			standalone computation from a single set of observations, via
 *  		e.g.: {@link #computeLocal(int[])},
 *  		{@link #computeAverageLocal(int[])} etc.</li>
 * 		<li>
 * 		Return to step 2 to re-use the calculator on a new data set.
 * 		</li>
 * 	</ol>
 * 
 * <p>TODO Inherit from {@link SingleAgentMeasureDiscreteInContextOfPastCalculator}
 *  as {@link ActiveInformationCalculatorDiscrete} does; Tidy up the Javadocs for
 *  the methods, which are somewhat preliminary</p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Bialek, W., Nemenman, I., and Tishby, N.,
 *  <a href="http://dx.doi.org/10.1016/S0378-4371(01)00444-7">
 * 	"Complexity through nonextensivity"</a>,
 *  Physica A, 302, 89-99. (2001).</li>
 * 	<li>J. P. Crutchfield, D. P. Feldman,
 *  <a href="http://dx.doi.org/10.1063/1.1530990">
 * 	"Regularities Unseen, Randomness Observed: Levels of Entropy Convergence"</a>,
 *  Chaos, Vol. 13, No. 1. (2003), pp. 25-54.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class PredictiveInformationCalculatorDiscrete {

	private double average = 0.0;
	private double max = 0.0;
	private double min = 0.0;
	private int observations = 0;
	private int k = 0; // history/future block length k.
	private int numDiscreteValues = 0; // number of individual states.
	private int[][]	jointCount = null; // Count for (x[t+1], x[t]) tuples
	private int[] prevCount = null; // Count for x[t]
	private int[] nextCount = null; // Count for x[t+1]
	private int[] maxShiftedValue = null; // states * (base^(history-1))
	
	private int base_power_k = 0;
	private double log_base = 0;
	
	/**
	 * User was formerly forced to create new instances through this factory method.
	 * Retained for backwards compatibility.
	 * 
	 * @param numDiscreteValues Number of discrete values (e.g. 2 for binary states)
	 * @param blockLength
	 * @deprecated
	 * @return
	 */
	public static PredictiveInformationCalculatorDiscrete newInstance(int numDiscreteValues, int blockLength) {
		return new PredictiveInformationCalculatorDiscrete(numDiscreteValues, blockLength);
	}
	
	/**
	 * Construct a new instance
	 * 
	 * @param numDiscreteValues number of quantisation levels for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param blockLength embedded history length of the past and future to use -
	 *        this is k in Schreiber's notation.
	 */
	public PredictiveInformationCalculatorDiscrete(int numDiscreteValues, int blockLength) {
		super();
		
		this.numDiscreteValues = numDiscreteValues;
		k = blockLength;
		base_power_k = MathsUtils.power(numDiscreteValues, k);
		log_base = Math.log(numDiscreteValues);
		
		if (blockLength < 1) {
			throw new RuntimeException("History k " + blockLength + " is not >= 1 for Predictive information calculator");
		}
		if (k > Math.log(Integer.MAX_VALUE) / log_base) {
			throw new RuntimeException("Base and history combination too large");
		}

		// Create storage for counts of observations
		jointCount = new int[base_power_k][base_power_k];
		prevCount = new int[base_power_k];
		nextCount = new int[base_power_k];

		// Create constants for tracking prevValues and nextValues
		maxShiftedValue = new int[numDiscreteValues];
		for (int v = 0; v < numDiscreteValues; v++) {
			maxShiftedValue[v] = v * MathsUtils.power(numDiscreteValues, k-1);
		}
	}

	/**
	 * Initialise calculator, preparing to take observation sets in
	 * Should be called prior to any of the addObservations() methods.
	 * You can reinitialise without needing to create a new object.
	 *
	 */
	public void initialise(){
		average = 0.0;
		max = 0.0;
		min = 0.0;
		observations = 0;
		
		MatrixUtils.fill(jointCount, 0);
		MatrixUtils.fill(prevCount, 0);
		MatrixUtils.fill(nextCount, 0);
	}
		
	/**
 	 * Add observations in to our estimates of the pdfs.
	 *
	 * @param timeSeries time series of agent states
	 */
	public void addObservations(int timeSeries[]) {
		int timeSteps = timeSeries.length;
		// increment the count of observations:
		// (we miss out on k observations at the start and k-1 at the end)
		if (timeSteps - k - (k-1) <= 0) {
			// Nothing to do
			return;
		}
		observations += (timeSteps - k - (k-1)); 
		
		// Initialise and store the current previous value 
		//  and next values (next val set for t = k-1)
		int prevVal = 0;
		int nextVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= numDiscreteValues;
			prevVal += timeSeries[p];
			nextVal *= numDiscreteValues;
			nextVal += timeSeries[k-1+p];
		}
		
		// 1. Count the tuples observed
		for (int t = k; t < timeSteps - (k-1); t++) {
			// Update the next value:
			nextVal -= maxShiftedValue[timeSeries[t-1]];
			nextVal *= numDiscreteValues;
			nextVal += timeSeries[k-1+t];
			// Update the counts
			jointCount[nextVal][prevVal]++;
			prevCount[prevVal]++;
			nextCount[nextVal]++;
			// Update the previous value:
			prevVal -= maxShiftedValue[timeSeries[t-k]];
			prevVal *= numDiscreteValues;
			prevVal += timeSeries[t];
		}		
	}

	/**
 	 * Add observations in to our estimates of the pdfs.
 	 * This call suitable only for homogeneous agents, as all
 	 *  agents will contribute to single pdfs.
	 *
	 * @param states 1st index is time, 2nd index is agent number
	 */
	public void addObservations(int states[][]) {
		int rows = states.length;
		int columns = states[0].length;
		// (we miss out on k observations at the start and k-1 at the end)
		if (rows - k - (k-1) <= 0) {
			// Nothing to do
			return;
		}
		observations += (rows - k - (k-1))*columns; 
		
		// Initialise and store the current previous and
		//  next val (next val set for t = k-1) value for each column
		int[] prevVal = new int[columns]; 
		int[] nextVal = new int[columns]; 
		for (int c = 0; c < columns; c++) {
			prevVal[c] = 0;
			nextVal[c] = 0;
			for (int p = 0; p < k; p++) {
				prevVal[c] *= numDiscreteValues;
				prevVal[c] += states[p][c];
				nextVal[c] *= numDiscreteValues;
				nextVal[c] += states[k-1+p][c];
			}
		}
		
		// 1. Count the tuples observed
		for (int r = k; r < rows - (k-1); r++) {
			for (int c = 0; c < columns; c++) {
				// Update the next value:
				nextVal[c] -= maxShiftedValue[states[r-1][c]];
				nextVal[c] *= numDiscreteValues;
				nextVal[c] += states[k-1+r][c];
				// Update the counts
				jointCount[nextVal[c]][prevVal[c]]++;
				prevCount[prevVal[c]]++;
				nextCount[nextVal[c]]++;
				// Update the previous value:
				prevVal[c] -= maxShiftedValue[states[r-k][c]];
				prevVal[c] *= numDiscreteValues;
				prevVal[c] += states[r][c];
			}
		}		
	}

	/**
 	 * Add observations in to our estimates of the pdfs.
 	 * This call suitable only for homogeneous agents, as all
 	 *  agents will contribute to single pdfs.
	 *
	 * @param states 1st index is time, 2nd and 3rd index give the 2D agent number
	 */
	public void addObservations(int states[][][]) {
		int timeSteps = states.length;
		// (we miss out on k observations at the start and k-1 at the end)
		if (timeSteps - k - (k-1) <= 0) {
			// Nothing to do
			return;
		}
		int agentRows = states[0].length;
		if (agentRows == 0) {
			return;
		}
		int agentColumns = states[0][0].length;
		// increment the count of observations:
		observations += (timeSteps - k - (k-1)) * agentRows * agentColumns; 

		// Initialise and store the current previous and
		// next (next val set for t = k-1) value for each column
		int[][] prevVal = new int[agentRows][agentColumns];
		int[][] nextVal = new int[agentRows][agentColumns];
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
				prevVal[r][c] = 0;
				nextVal[r][c] = 0;
				for (int p = 0; p < k; p++) {
					prevVal[r][c] *= numDiscreteValues;
					prevVal[r][c] += states[p][r][c];
					nextVal[r][c] *= numDiscreteValues;
					nextVal[r][c] += states[k-1+p][r][c];
				}
			}
		}
		
		// 1. Count the tuples observed
		for (int t = k; t < timeSteps - (k-1); t++) {
			for (int r = 0; r < agentRows; r++) {
				for (int c = 0; c < agentColumns; c++) {
					// Update the next value:
					nextVal[r][c] -= maxShiftedValue[states[t-1][r][c]];
					nextVal[r][c] *= numDiscreteValues;
					nextVal[r][c] += states[k-1+t][r][c];
					// Update the counts
					jointCount[nextVal[r][c]][prevVal[r][c]]++;
					prevCount[prevVal[r][c]]++;
					nextCount[nextVal[r][c]]++;
					// Update the previous value:
					prevVal[r][c] -= maxShiftedValue[states[t-k][r][c]];
					prevVal[r][c] *= numDiscreteValues;
					prevVal[r][c] += states[t][r][c];
				}
			}
		}		
	}

	/**
 	 * Add observations for a single agent of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to addObservations(int states[][])
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states 1st index is time, 2nd index is agent number
	 */
	public void addObservations(int states[][], int col) {
		int rows = states.length;
		// (we miss out on k observations at the start and k-1 at the end)
		if (rows - k - (k-1) <= 0) {
			// Nothing to do
			return;
		}
		observations += (rows - k - (k-1)); 
		
		// Initialise and store the current previous and
		//  next value (next val set for t = k-1) for each column
		int prevVal = 0; 
		int nextVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= numDiscreteValues;
			prevVal += states[p][col];
			nextVal *= numDiscreteValues;
			nextVal += states[k-1+p][col];
		}
		
		// 1. Count the tuples observed
		for (int r = k; r < rows - (k-1); r++) {
			// Update the next value:
			nextVal -= maxShiftedValue[states[r-1][col]];
			nextVal *= numDiscreteValues;
			nextVal += states[k-1+r][col];
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			jointCount[nextVal][prevVal]++;
			prevCount[prevVal]++;
			nextCount[nextVal]++;					
			// Update the previous value:
			prevVal -= maxShiftedValue[states[r-k][col]];
			prevVal *= numDiscreteValues;
			prevVal += states[r][col];
		}
	}

	/**
 	 * Add observations for a single agent of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to addObservations(int states[][])
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states 1st index is time, 2nd and 3rd index give the 2D agent number
	 */
	public void addObservations(int states[][][], int agentIndex1, int agentIndex2) {
		int timeSteps = states.length;
		// (we miss out on k observations at the start and k-1 at the end)
		if (timeSteps - k - (k-1) <= 0) {
			// Nothing to do
			return;
		}
		// increment the count of observations:
		observations += (timeSteps - k - (k-1)); 
		
		// Initialise and store the current previous and
		//  next value (next val set for t = k-1) for each column
		int prevVal = 0; 
		int nextVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= numDiscreteValues;
			prevVal += states[p][agentIndex1][agentIndex2];
			nextVal *= numDiscreteValues;
			nextVal += states[k-1+p][agentIndex1][agentIndex2];
		}
		
		// 1. Count the tuples observed
		for (int t = k; t < timeSteps - (k-1); t++) {
			// Update the next value:
			nextVal -= maxShiftedValue[states[t-1][agentIndex1][agentIndex2]];
			nextVal *= numDiscreteValues;
			nextVal += states[k-1+t][agentIndex1][agentIndex2];
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			jointCount[nextVal][prevVal]++;
			prevCount[prevVal]++;
			nextCount[nextVal]++;
			// Update the previous value:
			prevVal -= maxShiftedValue[states[t-k][agentIndex1][agentIndex2]];
			prevVal *= numDiscreteValues;
			prevVal += states[t][agentIndex1][agentIndex2];
		}
	}

	/**
	 * Returns the average predictive information from
	 *  the observed values which have been passed in previously. 
	 *  
	 * @return
	 */
	public double computeAverageLocalOfObservations() {
		double mi = 0.0;
		double miCont = 0.0;

		max = 0;
		min = 0;
		for (int nextVal = 0; nextVal < base_power_k; nextVal++) {
			// compute p_next
			double p_next = (double) nextCount[nextVal] / (double) observations;
			for (int prevVal = 0; prevVal < base_power_k; prevVal++) {
				// compute p_prev
				double p_prev = (double) prevCount[prevVal] / (double) observations;
				// compute p(prev, next)
				double p_joint = (double) jointCount[nextVal][prevVal] / (double) observations;
				// Compute MI contribution:
				if (p_joint > 0.0) {
					double logTerm = p_joint / (p_next * p_prev);
					double localValue = Math.log(logTerm) / log_base;
					miCont = p_joint * localValue;
					if (localValue > max) {
						max = localValue;
					} else if (localValue < min) {
						min = localValue;
					}
				} else {
					miCont = 0.0;
				}
				mi += miCont;
			}
		}
		
		average = mi;
		return mi;
	}
	
	/**
	 * Computes local predictive info for the given values
	 *  
	 * @param next joint state of future (see code for {@link #addObservations(int[])}
	 *   for how the joint integer value is computed)
	 * @param past joint state of past (see code for {@link #addObservations(int[])}
	 *   for how the joint integer value is computed)
	 * @return
	 */
	public double computeLocalFromPreviousObservations(int next, int past){
		double logTerm = ( (double) jointCount[next][past] ) /
		  ( (double) nextCount[next] *
			(double) prevCount[past] );
		logTerm *= (double) observations;
		return Math.log(logTerm) / log_base;
	}

	/**
	 * Computes local predictive information for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 *  
	 * @param timeSeries time series of values
	 * @return array of local predictive information values
	 */
	public double[] computeLocalFromPreviousObservations(int timeSeries[]){
		int timeSteps = timeSeries.length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localPredictive = new double[timeSteps];
		if (timeSteps < k + (k-1)) {
			// Nothing to do
			return localPredictive;
		}

		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous and
		//  next value (next val set for t = k-1)
		int prevVal = 0;
		int nextVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= numDiscreteValues;
			prevVal += timeSeries[p];
			nextVal *= numDiscreteValues;
			nextVal += timeSeries[k-1+p];
		}

		double logTerm = 0.0;
		for (int t = k; t < timeSteps - (k-1); t++) {
			// Update the next value:
			nextVal -= maxShiftedValue[timeSeries[t-1]];
			nextVal *= numDiscreteValues;
			nextVal += timeSeries[k-1+t];
			logTerm = ( (double) jointCount[nextVal][prevVal] ) /
			  		  ( (double) nextCount[nextVal] *
			  			(double) prevCount[prevVal] );
			// Now account for the fact that we've
			//  just used counts rather than probabilities,
			//  and we've got two counts on the bottom
			//  but one count on the top:
			logTerm *= (double) observations;
			localPredictive[t] = Math.log(logTerm) / log_base;
			average += localPredictive[t];
			if (localPredictive[t] > max) {
				max = localPredictive[t];
			} else if (localPredictive[t] < min) {
				min = localPredictive[t];
			}
			// Update the previous value:
			prevVal -= maxShiftedValue[timeSeries[t-k]];
			prevVal *= numDiscreteValues;
			prevVal += timeSeries[t];
		}
		average = average/(double) (timeSteps - k - (k-1));
		
		return localPredictive;
		
	}

	/**
	 * Computes local predictive information for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method to be used for homogeneous agents only, since
	 *  all observations are pooled together.
	 *  
	 * @param timeSeries 1st index is time, 2nd index is agent number
	 * @return
	 */
	public double[][] computeLocalFromPreviousObservations(int timeSeries[][]){
		int rows = timeSeries.length;
		int columns = timeSeries[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][] localPredictive = new double[rows][columns];
		if (rows < k + (k-1)) {
			// Nothing to do
			return localPredictive;
		}

		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous and
		//  next value (next val set for t = k-1) for each column
		int[] prevVal = new int[columns];
		int[] nextVal = new int[columns];
		for (int c = 0; c < columns; c++) {
			prevVal[c] = 0;
			nextVal[c] = 0;
			for (int p = 0; p < k; p++) {
				prevVal[c] *= numDiscreteValues;
				prevVal[c] += timeSeries[p][c];
				nextVal[c] *= numDiscreteValues;
				nextVal[c] += timeSeries[k-1+p][c];
			}
		}
		double logTerm = 0.0;
		for (int r = k; r < rows - (k-1); r++) {
			for (int c = 0; c < columns; c++) {
				// Update the next value:
				nextVal[c] -= maxShiftedValue[timeSeries[r-1][c]];
				nextVal[c] *= numDiscreteValues;
				nextVal[c] += timeSeries[k-1+r][c];
				logTerm = ( (double) jointCount[nextVal[c]][prevVal[c]] ) /
				  		  ( (double) nextCount[nextVal[c]] *
				  			(double) prevCount[prevVal[c]] );
				// Now account for the fact that we've
				//  just used counts rather than probabilities,
				//  and we've got two counts on the bottom
				//  but one count on the top:
				logTerm *= (double) observations;
				localPredictive[r][c] = Math.log(logTerm) / log_base;
				average += localPredictive[r][c];
				if (localPredictive[r][c] > max) {
					max = localPredictive[r][c];
				} else if (localPredictive[r][c] < min) {
					min = localPredictive[r][c];
				}
				// Update the previous value:
				prevVal[c] -= maxShiftedValue[timeSeries[r-k][c]];
				prevVal[c] *= numDiscreteValues;
				prevVal[c] += timeSeries[r][c];
			}
		}
		average = average/(double) (columns * (rows - k - (k-1)));
		
		return localPredictive;
		
	}
	
	/**
	 * Computes local active information storage for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method to be used for homogeneous agents only, since
	 *  all observations are pooled together.
	 *  
	 * @param timeSeries 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @return
	 */
	public double[][][] computeLocalFromPreviousObservations(int timeSeries[][][]){
		int timeSteps = timeSeries.length;
		int agentRows = timeSeries[0].length;
		int agentColumns = timeSeries[0][0].length;

		// Allocate for all rows even though we'll leave the first and
		// last ones as zeros
		double[][][] localPredictive = new double[timeSteps][agentRows][agentColumns];
		if (timeSteps < k + (k-1)) {
			// Nothing to do
			return localPredictive;
		}
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous and
		//  next value (next val set for t = k-1) for each agent
		int[][] prevVal = new int[agentRows][agentColumns]; 
		int[][] nextVal = new int[agentRows][agentColumns]; 
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
				prevVal[r][c] = 0;
				nextVal[r][c] = 0;
				for (int p = 0; p < k; p++) {
					prevVal[r][c] *= numDiscreteValues;
					prevVal[r][c] += timeSeries[p][r][c];
					nextVal[r][c] *= numDiscreteValues;
					nextVal[r][c] += timeSeries[k-1+p][r][c];
				}
			}
		}
		double logTerm = 0.0;
		for (int t = k; t < timeSteps - (k-1); t++) {
			for (int r = 0; r < agentRows; r++) {
				for (int c = 0; c < agentColumns; c++) {
					// Update the next value:
					nextVal[r][c] -= maxShiftedValue[timeSeries[t-1][r][c]];
					nextVal[r][c] *= numDiscreteValues;
					nextVal[r][c] += timeSeries[k-1+t][r][c];
					logTerm = ( (double) jointCount[nextVal[r][c]][prevVal[r][c]] ) /
					  		  ( (double) nextCount[nextVal[r][c]] *
					  			(double) prevCount[prevVal[r][c]] );
					// Now account for the fact that we've
					//  just used counts rather than probabilities,
					//  and we've got two counts on the bottom
					//  but one count on the top:
					logTerm *= (double) observations;
					localPredictive[t][r][c] = Math.log(logTerm) / log_base;
					average += localPredictive[t][r][c];
					if (localPredictive[t][r][c] > max) {
						max = localPredictive[t][r][c];
					} else if (localPredictive[t][r][c] < min) {
						min = localPredictive[t][r][c];
					}
					// Update the previous value:
					prevVal[r][c] -= maxShiftedValue[timeSeries[t-k][r][c]];
					prevVal[r][c] *= numDiscreteValues;
					prevVal[r][c] += timeSeries[t][r][c];
				}
			}
		}
		average = average/(double) (agentRows * agentColumns * (timeSteps - k - (k-1)));
		
		return localPredictive;
	}

	/**
	 * Computes local predictive information for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method is suitable for heterogeneous agents since the 
	 *  user specifies which agent to take the observations from.
	 *  
	 * @param states 1st index is time, 2nd index is agent number
	 * @param col gives the column index identifying the agent
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations(int states[][], int col){
		int rows = states.length;
		//int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localPredictive = new double[rows];
		if (rows < k + (k-1)) {
			// Nothing to do
			return localPredictive;
		}
		average = 0;
		max = 0;
		min = 0;
		
		// Initialise and store the current previous and 
		//  next value (next val set for t = k-1) for each column
		int prevVal = 0; 
		int nextVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= numDiscreteValues;
			prevVal += states[p][col];
			nextVal *= numDiscreteValues;
			nextVal += states[k-1+p][col];
		}
		double logTerm = 0.0;
		for (int r = k; r < rows - (k-1); r++) {
			// Update the next value:
			nextVal -= maxShiftedValue[states[r-1][col]];
			nextVal *= numDiscreteValues;
			nextVal += states[k-1+r][col];
			logTerm = ( (double) jointCount[nextVal][prevVal] ) /
			  		  ( (double) nextCount[nextVal] *
			  			(double) prevCount[prevVal] );
			// Now account for the fact that we've
			//  just used counts rather than probabilities,
			//  and we've got two counts on the bottom
			//  but one count on the top:
			logTerm *= (double) observations;
			localPredictive[r] = Math.log(logTerm) / log_base;
			average += localPredictive[r];
			if (localPredictive[r] > max) {
				max = localPredictive[r];
			} else if (localPredictive[r] < min) {
				min = localPredictive[r];
			}
			// Update the previous value:
			prevVal -= maxShiftedValue[states[r-k][col]];
			prevVal *= numDiscreteValues;
			prevVal += states[r][col];
		}
		average = average/(double) (rows - k - (k-1));
		
		return localPredictive;
		
	}
	
	/**
	 * Computes local predictive information for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method is suitable for heterogeneous agents, since the
	 *  relevant agent is identified
	 *  
	 * @param timeSeries 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @param agentIndex1 gives the first index identifying the agent
	 * @param agentIndex2 gives the second index identifying the agent
	 * @return
	 */
	public double[] computeLocalFromPreviousObservations(int timeSeries[][][], int agentIndex1, int agentIndex2){
		int timeSteps = timeSeries.length;
		//int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localPredictive = new double[timeSteps];
		if (timeSteps < k + (k-1)) {
			// Nothing to do
			return localPredictive;
		}
		average = 0;
		max = 0;
		min = 0;
		
		// Initialise and store the current previous and 
		//  next value (next val set for t = k-1) for each column
		int prevVal = 0;
		int nextVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= numDiscreteValues;
			prevVal += timeSeries[p][agentIndex1][agentIndex2];
			nextVal *= numDiscreteValues;
			nextVal += timeSeries[k-1+p][agentIndex1][agentIndex2];
		}
		double logTerm = 0.0;
		for (int t = k; t < timeSteps - (k-1); t++) {
			// Update the next value:
			nextVal -= maxShiftedValue[timeSeries[t-1][agentIndex1][agentIndex2]];
			nextVal *= numDiscreteValues;
			nextVal += timeSeries[k-1+t][agentIndex1][agentIndex2];
			logTerm = ( (double) jointCount[nextVal][prevVal] ) /
			  		  ( (double) nextCount[nextVal] *
			  			(double) prevCount[prevVal] );
			// Now account for the fact that we've
			//  just used counts rather than probabilities,
			//  and we've got two counts on the bottom
			//  but one count on the top:
			logTerm *= (double) observations;
			localPredictive[t] = Math.log(logTerm) / log_base;
			average += localPredictive[t];
			if (localPredictive[t] > max) {
				max = localPredictive[t];
			} else if (localPredictive[t] < min) {
				min = localPredictive[t];
			}
			// Update the previous value:
			prevVal -= maxShiftedValue[timeSeries[t-k][agentIndex1][agentIndex2]];
			prevVal *= numDiscreteValues;
			prevVal += timeSeries[t][agentIndex1][agentIndex2];
		}
		average = average/(double) (timeSteps - k - (k-1));
		
		return localPredictive;
		
	}

	/**
	 * Standalone routine to 
	 * compute local predictive information from a time series.
	 * Return an array of local values.
	 * First k and last k-1 rows are zeros.
	 * 
	 * @param timeSeries time series array of states
	 * @return time series of local predictive information values
	 */
	public double[] computeLocal(int timeSeries[]) {
		
		initialise();
		addObservations(timeSeries);
		return computeLocalFromPreviousObservations(timeSeries);
	}

	/**
	 * Standalone routine to 
	 * compute local predictive information across a 2D spatiotemporal
	 *  time series of the observations of homogeneous agents.
	 * This method to be called for homogeneous agents only, as all 
	 *  observations are pooled together.
	 * First k and last k-1 rows are zeros.
	 * 
	 * @param timeSeries 2D array of states (first index is time, second is variable id)
	 * @return multivariate time series of local predictive information values
	 *  (1st index is time, second is variable id)
	 */
	public double[][] computeLocal(int timeSeries[][]) {
		
		initialise();
		addObservations(timeSeries);
		return computeLocalFromPreviousObservations(timeSeries);
	}
	
	/**
	 * Standalone routine to 
	 * compute local predictive information across a 3D spatiotemporal
	 *  array of the states of homogeneous agents
	 * This method to be called for homogeneous agents only, as all 
	 *  observations are pooled together.
	 * First k and last k-1 rows are zeros.
	 * 
	 * @param timeSeries 2D array of states (first index is time, second
	 *  and third are variable ids)
	 * @return multivariate time series of local predictive information values
	 *  (1st index is time, second and third are variable ids)
	 */
	public double[][][] computeLocal(int timeSeries[][][]) {
		
		initialise();
		addObservations(timeSeries);
		return computeLocalFromPreviousObservations(timeSeries);
	}

	/**
	 * Standalone routine to 
	 * compute average predictive information across a
	 *  time series.
	 * 
	 * @param timeSeries time series array of states
	 * @return average predictive information values
	 */
	public double computeAverageLocal(int timeSeries[]) {
		
		initialise();
		addObservations(timeSeries);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute average predictive information across a 2D spatiotemporal
	 *  array of the states of homogeneous agents.
	 * This method to be called for homogeneous agents only, as all 
	 *  observations are pooled together.
	 * 
	 * @param timeSeries 2D array of states (first index is time, second is variable id)
	 * @return average predictive information value
	 */
	public double computeAverageLocal(int timeSeries[][]) {
		
		initialise();
		addObservations(timeSeries);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute average predictive information across a 3D spatiotemporal
	 *  array of the states of homogeneous agents
	 * This method to be called for homogeneous agents only, as all 
	 *  observations are pooled together.
	 * 
	 * @param timeSeries 2D array of states (first index is time, second
	 *  and third are variable ids)
	 * @return average predictive information value
	 */
	public double computeAverageLocal(int timeSeries[][][]) {
		
		initialise();
		addObservations(timeSeries);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute local predictive information for one agent in a 2D spatiotemporal
	 *  array of the states of agents
	 * First k and last k-1 rows are zeros.
	 * 
	 * This method should be used for heterogeneous agents.
	 * 
	 * @param timeSeries 2D array of states (first index is time, second is variable id)
	 * @param col index of the agent whose local values should be calculated
	 * @return time series of local predictive information values for the given agent id
	 */
	public double[] computeLocal(int timeSeries[][], int col) {
		
		initialise();
		addObservations(timeSeries, col);
		return computeLocalFromPreviousObservations(timeSeries, col);
	}
	
	/**
	 * Standalone routine to 
	 * compute local predictive information for one agent in a 3D spatiotemporal
	 *  array of the states of agents
	 * First k and last k-1 rows are zeros.
	 * This method should be used for heterogeneous agents.
	 * 
	 * @param timeSeries 3D array of states (first index is time, second 
	 *  and third are variable id)
	 * @param agentIndex1 row index of agent
	 * @param agentIndex2 column index of agent
	 * @return time series of local predictive information values for the given agent id
	 */
	public double[] computeLocal(int timeSeries[][][], int agentIndex1, int agentIndex2) {
		
		initialise();
		addObservations(timeSeries, agentIndex1, agentIndex2);
		return computeLocalFromPreviousObservations(timeSeries, agentIndex1, agentIndex2);
	}

	/**
	 * Standalone routine to 
	 * compute average local predictive information 
	 * for a single agent in a 2D spatiotemporal time series.
	 * This method suitable for heterogeneous agents, since the agent
	 *  to be investigated is named.
	 * 
	 * @param timeSeries 2D array of states (first index is time, second is variable id)
	 * @param col index of the agent whose local values should be calculated
	 * @return average predictive info for the given agent id
	 */
	public double computeAverageLocal(int timeSeries[][], int col) {
		
		initialise();
		addObservations(timeSeries, col);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute average local active information storage 
	 * for a single agent in a 3D spatiotemporal time series.
	 * This method suitable for heterogeneous agents, since the agent
	 *  to be investigated is named.
	 * 
	 * @param timeSeries 3D array of states (first index is time, second 
	 *  and third are variable id)
	 * @param agentIndex1 row index of agent
	 * @param agentIndex2 column index of agent
	 * @return average predictive info for the given agent id
	 */
	public double computeAverageLocal(int timeSeries[][][], int agentIndex1, int agentIndex2) {
		
		initialise();
		addObservations(timeSeries, agentIndex1, agentIndex2);
		return computeAverageLocalOfObservations();
	}

	/**
	 * 
	 * @return the last predictive information computed since the 
	 *  previous {@link #initialise()} call.
	 */
	public double getLastAverage() {
		return average;
	}

	/**
	 * 
	 * @return the last maximum local predictive information computed since the 
	 *  previous {@link #initialise()} call.
	 */
	public double getLastMax() {
		return max;
	}

	/**
	 * 
	 * @return the last minimum local predictive information computed since the 
	 *  previous {@link #initialise()} call.
	 */
	public double getLastMin() {
		return min;
	}
	
	/**
	 * Writes the current probability distribution functions 
	 *  
	 * @return
	 */
	public void writePdfs() {
		double mi = 0.0;
		double miCont = 0.0;

		System.out.println("nextVal p(next) prevVal p(prev) p(joint) logTerm localVal");
		for (int nextVal = 0; nextVal < base_power_k; nextVal++) {
			// compute p_next
			double p_next = (double) nextCount[nextVal] / (double) observations;
			for (int prevVal = 0; prevVal < base_power_k; prevVal++) {
				// compute p_prev
				double p_prev = (double) prevCount[prevVal] / (double) observations;
				// compute p(prev, next)
				double p_joint = (double) jointCount[nextVal][prevVal] / (double) observations;
				// Compute MI contribution:
				if (p_joint > 0.0) {
					double logTerm = p_joint / (p_next * p_prev);
					double localValue = Math.log(logTerm) / log_base;
					miCont = p_joint * localValue;
					System.out.println(String.format("%7d    %.2f %7d    %.2f     %.2f    %.2f     %.2f",
							nextVal, p_next, prevVal, p_prev, p_joint, logTerm, localValue));
				} else {
					miCont = 0.0;
					System.out.println(String.format("%7d    %.2f %7d    %.2f     %.2f    %.2f     %.2f",
							nextVal, p_next, prevVal, p_prev, p_joint, 0.0, 0.0));
				}
				mi += miCont;
			}
		}
		System.out.println("Average is " + mi);
		
		return;
	}

	/**
	 * Utility function to compute the combined past values of x up to and including time step t
	 *  (i.e. (x_{t-k+1}, ... ,x_{t-1},x_{t}))
	 * 
	 * @param x time series
	 * @param t time step to compute joint past value up to
	 * @return joint discrete value computed from the (x_{t-k+1}, ... ,x_{t-1},x_{t})
	 */
	public int computePastValue(int[] x, int t) {
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= numDiscreteValues;
			pastVal += x[t - k + 1 + p];
		}
		return pastVal;
	}

	/**
	 * Utility function to compute the combined past values of x up to
	 * and including time step t for a given variable number within x
	 *  (i.e. (x_{t-k+1,i}, ... ,x_{t-1,i},x_{t,i}))
	 * 
	 * @param x 2D multivariate time series
	 * @param i column or agent number within x.
	 * @param t time step to compute joint past value up to
	 * @return joint discrete value computed from the (x_{t-k+1,i}, ... ,x_{t-1,i},x_{t,i}))
	 */
	public int computePastValue(int[][] x, int i, int t) {
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= numDiscreteValues;
			pastVal += x[t - k + 1 + p][i];
		}
		return pastVal;
	}

	/**
	 * Utility function to compute the combined past values of x up to
	 * and including time step t for a given variable number within x
	 *  (i.e. (x_{t-k+1,i,j}, ... ,x_{t-1,i,j},x_{t,i,j}))
	 * 
	 * @param x 3D multivariate time series
	 * @param i 1st index to identify agent within x.
	 * @param i 2nd index to identify agent within x.
	 * @param t time step to compute joint past value up to
	 * @return joint discrete value computed from the (x_{t-k+1,i}, ... ,x_{t-1,i},x_{t,i}))
	 */
	public int computePastValue(int[][][] x, int i, int j, int t) {
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= numDiscreteValues;
			pastVal += x[t - k + 1 + p][i][j];
		}
		return pastVal;
	}
}
