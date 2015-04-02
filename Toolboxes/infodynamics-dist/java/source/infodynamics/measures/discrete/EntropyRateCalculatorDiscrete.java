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
 * <p>Entropy rate calculator for univariate discrete (int[]) data
 * (ie computes entropy over blocks of consecutive states in time).
 * Implements entropy rate as entropy of next state
 * conditional on the embedded past (as per the alternative
 * definition used by Crutchfield and Feldman, see below)
 * rather than the limiting rate of block entropy over block size.</p>
 *
 * <p>Usage of the class is intended to follow this paradigm:</p>
 * <ol>
 * 		<li>Construct the calculator: {@link #EntropyRateCalculatorDiscrete(int, int)};</li>
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
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. M. Cover and J. A. Thomas, 'Elements of Information
Theory' (John Wiley & Sons, New York, 1991).</li>
 * 	<li>J. P. Crutchfield, D. P. Feldman,
 *  <a href="http://dx.doi.org/10.1063/1.1530990">
 * 	"Regularities Unseen, Randomness Observed: Levels of Entropy Convergence"</a>,
 *  Chaos, Vol. 13, No. 1. (2003), pp. 25-54.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class EntropyRateCalculatorDiscrete extends SingleAgentMeasureDiscreteInContextOfPastCalculator {

	/**
	 * User was formerly forced to create new instances through this factory method.
	 * Retained for backwards compatibility.
	 * 
	 * @param base
	 * @param history
	 * @deprecated
	 * @return
	 */
	public static EntropyRateCalculatorDiscrete newInstance(int base, int history) {
		return new EntropyRateCalculatorDiscrete(base, history);
	}
	
	/**
	 * Construct a new instance
	 * 
	 * @param base number of symbols for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param history embedded history length of the destination to condition on -
	 *        this is k in Schreiber's notation.
	 */
	public EntropyRateCalculatorDiscrete(int base, int history) {
		super(base, history);
	}		
	
	@Override
	public void addObservations(int[] states) {
		int rows = states.length;
		// increment the count of observations:
		observations += (rows - k); 
		
		// Initialise and store the current previous value for each column
		int prevVal = 0; 
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p];
		}
		
		// 1. Count the tuples observed
		int nextVal;
		for (int r = k; r < rows; r++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			nextVal = states[r];
			nextPastCount[nextVal][prevVal]++;
			pastCount[prevVal]++;
			// Update the previous value:
			prevVal -= maxShiftedValue[states[r-k]];
			prevVal *= base;
			prevVal += states[r];
		}
	}

	@Override
	public void addObservations(int states[][]) {
		int rows = states.length;
		int columns = states[0].length;
		// increment the count of observations:
		observations += (rows - k)*columns; 
		
		// Initialise and store the current previous value for each column
		int[] prevVal = new int[columns]; 
		for (int c = 0; c < columns; c++) {
			prevVal[c] = 0;
			for (int p = 0; p < k; p++) {
				prevVal[c] *= base;
				prevVal[c] += states[p][c];
			}
		}
		
		// 1. Count the tuples observed
		int nextVal;
		for (int r = k; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				// Add to the count for this particular transition:
				// (cell's assigned as above)
				nextVal = states[r][c];
				nextPastCount[nextVal][prevVal[c]]++;
				pastCount[prevVal[c]]++;
				// Update the previous value:
				prevVal[c] -= maxShiftedValue[states[r-k][c]];
				prevVal[c] *= base;
				prevVal[c] += states[r][c];
			}
		}		
	}
	
	@Override
	public void addObservations(int states[][][]) {
		int timeSteps = states.length;
		if (timeSteps == 0) {
			return;
		}
		int agentRows = states[0].length;
		if (agentRows == 0) {
			return;
		}
		int agentColumns = states[0][0].length;
		// increment the count of observations:
		observations += (timeSteps - k) * agentRows * agentColumns; 
		
		// Initialise and store the current previous value for each column
		int[][] prevVal = new int[agentRows][agentColumns];
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
				prevVal[r][c] = 0;
				for (int p = 0; p < k; p++) {
					prevVal[r][c] *= base;
					prevVal[r][c] += states[p][r][c];
				}
			}
		}
		
		// 1. Count the tuples observed
		int nextVal;
		for (int t = k; t < timeSteps; t++) {
			for (int r = 0; r < agentRows; r++) {
				for (int c = 0; c < agentColumns; c++) {
					// Add to the count for this particular transition:
					// (cell's assigned as above)
					nextVal = states[t][r][c];
					nextPastCount[nextVal][prevVal[r][c]]++;
					pastCount[prevVal[r][c]]++;
					// Update the previous value:
					prevVal[r][c] -= maxShiftedValue[states[t-k][r][c]];
					prevVal[r][c] *= base;
					prevVal[r][c] += states[t][r][c];
				}
			}
		}		
	}

	@Override
	public void addObservations(int states[][], int col) {
		int rows = states.length;
		// increment the count of observations:
		observations += (rows - k); 
		
		// Initialise and store the current previous value for each column
		int prevVal = 0; 
		prevVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p][col];
		}
		
		// 1. Count the tuples observed
		int nextVal;
		for (int r = k; r < rows; r++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			nextVal = states[r][col];
			nextPastCount[nextVal][prevVal]++;
			pastCount[prevVal]++;
			// Update the previous value:
			prevVal -= maxShiftedValue[states[r-k][col]];
			prevVal *= base;
			prevVal += states[r][col];
		}
	}

	@Override
	public void addObservations(int states[][][], int agentIndex1, int agentIndex2) {
		int timeSteps = states.length;
		// increment the count of observations:
		observations += (timeSteps - k); 
		
		// Initialise and store the current previous value for this column
		int prevVal = 0; 
		prevVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p][agentIndex1][agentIndex2];
		}
		
		// 1. Count the tuples observed
		int nextVal;
		for (int r = k; r < timeSteps; r++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			nextVal = states[r][agentIndex1][agentIndex2];
			nextPastCount[nextVal][prevVal]++;
			pastCount[prevVal]++;
			// Update the previous value:
			prevVal -= maxShiftedValue[states[r-k][agentIndex1][agentIndex2]];
			prevVal *= base;
			prevVal += states[r][agentIndex1][agentIndex2];
		}
	}

	@Override
	public double computeAverageLocalOfObservations() {
		double entRate = 0.0;
		double entRateCont = 0.0;

		max = 0;
		min = 0;
		double logTerm = 0;
		for (int nextVal = 0; nextVal < base; nextVal++) {
			for (int prevVal = 0; prevVal < base_power_k; prevVal++) {
				// compute p_prev
				double p_prev = (double) pastCount[prevVal] / (double) observations;
				// compute p(prev, next)
				double p_joint = (double) nextPastCount[nextVal][prevVal] / (double) observations;
				// Compute entropy rate contribution:
				if (p_joint > 0.0) {
					logTerm = p_joint / p_prev;
					// Entropy rate takes the negative log:
					double localValue = - Math.log(logTerm) / log_2;
					entRateCont = p_joint * localValue;
					if (localValue > max) {
						max = localValue;
					} else if (localValue < min) {
						min = localValue;
					}
				} else {
					entRateCont = 0.0;
				}
				entRate += entRateCont;
			}
		}
		
		average = entRate;
		return entRate;
	}
	
	@Override
	public double[] computeLocalFromPreviousObservations(int[] states) {
		int rows = states.length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localEntRate = new double[rows];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int prevVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p];
		}
		
		int nextVal;
		double logTerm = 0.0;
		for (int r = k; r < rows; r++) {
			nextVal = states[r];
			logTerm = ( (double) nextPastCount[nextVal][prevVal] ) /
			  		  ( (double) pastCount[prevVal] );
			// Entropy rate takes the negative log:
			localEntRate[r] = - Math.log(logTerm) / log_2;
			average += localEntRate[r];
			if (localEntRate[r] > max) {
				max = localEntRate[r];
			} else if (localEntRate[r] < min) {
				min = localEntRate[r];
			}
			// Update the previous value:
			prevVal -= maxShiftedValue[states[r-k]];
			prevVal *= base;
			prevVal += states[r];
		}
		average = average/(double) (rows - k);
		
		return localEntRate;
	}

	@Override
	public double[][] computeLocalFromPreviousObservations(int states[][]){
		int rows = states.length;
		int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][] localEntRate = new double[rows][columns];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int[] prevVal = new int[columns]; 
		for (int c = 0; c < columns; c++) {
			prevVal[c] = 0;
			for (int p = 0; p < k; p++) {
				prevVal[c] *= base;
				prevVal[c] += states[p][c];
			}
		}
		int nextVal;
		double logTerm = 0.0;
		for (int r = k; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				nextVal = states[r][c];
				logTerm = ( (double) nextPastCount[nextVal][prevVal[c]] ) /
				  		  ( (double) pastCount[prevVal[c]] );
				// Entropy rate takes the negative log:
				localEntRate[r][c] = - Math.log(logTerm) / log_2;
				average += localEntRate[r][c];
				if (localEntRate[r][c] > max) {
					max = localEntRate[r][c];
				} else if (localEntRate[r][c] < min) {
					min = localEntRate[r][c];
				}
				// Update the previous value:
				prevVal[c] -= maxShiftedValue[states[r-k][c]];
				prevVal[c] *= base;
				prevVal[c] += states[r][c];
			}
		}
		average = average/(double) (columns * (rows - k));
		
		return localEntRate;		
	}
	
	@Override
	public double[][][] computeLocalFromPreviousObservations(int states[][][]){
		int timeSteps = states.length;
		int agentRows = states[0].length;
		int agentColumns = states[0][0].length;

		// Allocate for all time steps even though we'll leave the first ones as zeros
		double[][][] localEntRate = new double[timeSteps][agentRows][agentColumns];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int[][] prevVal = new int[agentRows][agentColumns];
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
				prevVal[r][c] = 0;
				for (int p = 0; p < k; p++) {
					prevVal[r][c] *= base;
					prevVal[r][c] += states[p][r][c];
				}
			}
		}
		
		int nextVal;
		double logTerm = 0.0;
		for (int t = k; t < timeSteps; t++) {
			for (int r = 0; r < agentRows; r++) {
				for (int c = 0; c < agentColumns; c++) {
					nextVal = states[t][r][c];
					logTerm = ( (double) nextPastCount[nextVal][prevVal[r][c]] ) /
					  		  ( (double) pastCount[prevVal[r][c]] );
					// Entropy rate takes the negative log:
					localEntRate[t][r][c] = - Math.log(logTerm) / log_2;
					average += localEntRate[t][r][c];
					if (localEntRate[t][r][c] > max) {
						max = localEntRate[t][r][c];
					} else if (localEntRate[t][r][c] < min) {
						min = localEntRate[t][r][c];
					}
					// Update the previous value:
					prevVal[r][c] -= maxShiftedValue[states[t-k][r][c]];
					prevVal[r][c] *= base;
					prevVal[r][c] += states[t][r][c];
				}
			}
		}
		average = average/(double) (agentRows * agentColumns * (timeSteps - k));
		
		return localEntRate;
	}

	@Override
	public double[] computeLocalFromPreviousObservations(int states[][], int col){
		int rows = states.length;
		//int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localEntRate = new double[rows];
		average = 0;
		max = 0;
		min = 0;
		
		// Initialise and store the current previous value for each column
		int prevVal = 0; 
		prevVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p][col];
		}
		int nextVal;
		double logTerm = 0.0;
		for (int r = k; r < rows; r++) {
			nextVal = states[r][col];
			logTerm = ( (double) nextPastCount[nextVal][prevVal] ) /
			  		  ( (double) pastCount[prevVal] );
			// Entropy rate takes the negative log:
			localEntRate[r] = - Math.log(logTerm) / log_2;
			average += localEntRate[r];
			if (localEntRate[r] > max) {
				max = localEntRate[r];
			} else if (localEntRate[r] < min) {
				min = localEntRate[r];
			}
			// Update the previous value:
			prevVal -= maxShiftedValue[states[r-k][col]];
			prevVal *= base;
			prevVal += states[r][col];
		}
		average = average/(double) (rows - k);
		
		return localEntRate;
		
	}

	@Override
	public double[] computeLocalFromPreviousObservations(int states[][][], int agentIndex1, int agentIndex2){
		int timeSteps = states.length;
		//int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localEntRate = new double[timeSteps];
		average = 0;
		max = 0;
		min = 0;
		
		// Initialise and store the current previous value for each column
		int prevVal = 0; 
		prevVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p][agentIndex1][agentIndex2];
		}
		int nextVal;
		double logTerm = 0.0;
		for (int t = k; t < timeSteps; t++) {
			nextVal = states[t][agentIndex1][agentIndex2];
			logTerm = ( (double) nextPastCount[nextVal][prevVal] ) /
			  		  ( (double) pastCount[prevVal] );
			// Entropy rate takes the negative log:
			localEntRate[t] = - Math.log(logTerm) / log_2;
			average += localEntRate[t];
			if (localEntRate[t] > max) {
				max = localEntRate[t];
			} else if (localEntRate[t] < min) {
				min = localEntRate[t];
			}
			// Update the previous value:
			prevVal -= maxShiftedValue[states[t-k][agentIndex1][agentIndex2]];
			prevVal *= base;
			prevVal += states[t][agentIndex1][agentIndex2];
		}
		average = average/(double) (timeSteps - k);
		
		return localEntRate;
		
	}
}
