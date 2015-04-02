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

import infodynamics.utils.MatrixUtils;

/**
 * <p>Entropy calculator for univariate discrete (int[]) data.</p>
 * 
 * <p>Usage of the class is intended to follow this paradigm:</p>
 * <ol>
 * 		<li>Construct the calculator: {@link #EntropyCalculatorDiscrete(int)};</li>
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
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class EntropyCalculatorDiscrete extends InfoMeasureCalculatorDiscrete
				implements SingleAgentMeasureDiscrete
{

	protected int[] stateCount = null; // Count for i[t]

	/**
	 * User was formerly forced to create new instances through this factory method.
	 * Retained for backwards compatibility.
	 * 
	 * @param base number of symbols for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param blocksize number of consecutive joint values to include
	 *  in the calculation.
	 * @deprecated
	 * @return a new EntropyCalculator
	 */
	public static EntropyCalculatorDiscrete newInstance(int base, int blocksize) {
		if (blocksize > 1) {
			return BlockEntropyCalculatorDiscrete.newInstance(blocksize, base);
		} else {
			return EntropyCalculatorDiscrete.newInstance(base);
		}
	}
	public static EntropyCalculatorDiscrete newInstance(int base) {
		return new EntropyCalculatorDiscrete(base);
	}
	
	/**
	 * Contruct a new instance
	 * 
	 * @param base number of quantisation levels for each variable.
	 *        E.g. binary variables are in base-2.
	 */
	public EntropyCalculatorDiscrete(int base) {

		super(base);
		
		// Create storage for counts of observations
		stateCount = new int[base];
	}
	
	@Override
	public void initialise(){
		super.initialise();
		MatrixUtils.fill(stateCount, 0);
	}
		
	@Override
	public void addObservations(int states[]) {
		int rows = states.length;
		// increment the count of observations:
		observations += rows; 
		
		// 1. Count the tuples observed
		for (int r = 0; r < rows; r++) {
			// Add to the count for this particular state:
			stateCount[states[r]]++;
		}		
	}

	@Override
	public void addObservations(int states[][]) {
		int rows = states.length;
		int columns = states[0].length;
		// increment the count of observations:
		observations += rows * columns; 
		
		// 1. Count the tuples observed
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				// Add to the count for this particular state:
				stateCount[states[r][c]]++;					
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
		observations += timeSteps * agentRows * agentColumns; 
		
		// 1. Count the tuples observed
		for (int t = 0; t < timeSteps; t++) {
			for (int i = 0; i < agentRows; i++) {
				for (int j = 0; j < agentColumns; j++) {
					// Add to the count for this particular state:
					stateCount[states[t][i][j]]++;
				}
			}
		}		
	}

	@Override
	public void addObservations(int states[][], int agentNumber) {
		int rows = states.length;
		// increment the count of observations:
		observations += rows; 
		
		// 1. Count the tuples observed
		for (int r = 0; r < rows; r++) {
			// Add to the count for this particular state:
			stateCount[states[r][agentNumber]]++;					
		}
	}

	@Override
	public void addObservations(int states[][][], int agentIndex1, int agentIndex2) {
		int timeSteps = states.length;
		// increment the count of observations:
		observations += timeSteps; 
		
		// 1. Count the tuples observed
		for (int r = 0; r < timeSteps; r++) {
			// Add to the count for this particular state:
			stateCount[states[r][agentIndex1][agentIndex2]]++;					
		}
	}

	/**
	 * Return the current count for the given value
	 * 
	 * @param stateVal given value
	 * @return count of observations of the given state
	 */
	public int getStateCount(int stateVal) {
		return stateCount[stateVal];
	}
	
	/**
	 * Return the current probability for the given value
	 * 
	 * @param stateVal given value
	 * @return probability of the given state
	 */
	public double getStateProbability(int stateVal) {
		return (double) stateCount[stateVal] / (double) observations;
	}

	@Override
	public double computeAverageLocalOfObservations() {
		double ent = 0.0;
		double entCont = 0.0;

		max = 0;
		min = 0;
		for (int stateVal = 0; stateVal < base; stateVal++) {
			// compute p_state
			double p_state = (double) stateCount[stateVal] / (double) observations;
			if (p_state > 0.0) {
				// Entropy takes the negative log:
				double localValue = - Math.log(p_state) / log_2;
				entCont = p_state * localValue;
				if (localValue > max) {
					max = localValue;
				} else if (localValue < min) {
					min = localValue;
				}
			} else {
				entCont = 0.0;
			}
			ent += entCont;
		}
		
		average = ent;
		return ent;
	}
	
	@Override
	public double[] computeLocalFromPreviousObservations(int states[]){
		int rows = states.length;
		
		double[] localEntropy = new double[rows];
		average = 0;
		max = 0;
		min = 0;
		for (int r = 0; r < rows; r++) {
			double p_state = (double) stateCount[states[r]] / (double) observations;
			// Entropy takes the negative log:
			localEntropy[r] = - Math.log(p_state) / log_2;
			average += localEntropy[r];
			if (localEntropy[r] > max) {
				max = localEntropy[r];
			} else if (localEntropy[r] < min) {
				min = localEntropy[r];
			}
		}
		average = average/(double) rows;
		
		return localEntropy;
		
	}
	
	@Override
	public double[][] computeLocalFromPreviousObservations(int states[][]){
		int rows = states.length;
		int columns = states[0].length;
		
		double[][] localEntropy = new double[rows][columns];
		average = 0;
		max = 0;
		min = 0;
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				double p_state = (double) stateCount[states[r][c]] / (double) observations;
				// Entropy takes the negative log:
				localEntropy[r][c] = - Math.log(p_state) / log_2;
				average += localEntropy[r][c];
				if (localEntropy[r][c] > max) {
					max = localEntropy[r][c];
				} else if (localEntropy[r][c] < min) {
					min = localEntropy[r][c];
				}
			}
		}
		average = average/(double) (columns * rows);
		
		return localEntropy;
		
	}
	
	@Override
	public double[][][] computeLocalFromPreviousObservations(int states[][][]){
		int timeSteps = states.length;
		int agentRows, agentColumns;
		if (timeSteps == 0) {
			agentRows = 0;
			agentColumns = 0;
		} else {
			agentRows = states[0].length;
			if (agentRows == 0) {
				agentColumns = 0;
			} else {
				agentColumns = states[0][0].length;
			}
		}

		double[][][] localEntropy = new double[timeSteps][agentRows][agentColumns];
		average = 0;
		max = 0;
		min = 0;
		for (int r = 0; r < timeSteps; r++) {
			for (int i = 0; i < agentRows; i++) {
				for (int j = 0; j < agentColumns; j++) {
					double p_state = (double) stateCount[states[r][i][j]] / (double) observations;
					// Entropy takes the negative log:
					localEntropy[r][i][j] = - Math.log(p_state) / log_2;
					average += localEntropy[r][i][j];
					if (localEntropy[r][i][j] > max) {
						max = localEntropy[r][i][j];
					} else if (localEntropy[r][i][j] < min) {
						min = localEntropy[r][i][j];
					}
				}
			}
		}
		average = average/(double) (agentRows * agentColumns * timeSteps);
		
		return localEntropy;
		
	}

	@Override
	public double[] computeLocalFromPreviousObservations(int states[][], int agentNumber){
		int rows = states.length;
		//int columns = states[0].length;

		double[] localEntropy = new double[rows];
		average = 0;
		max = 0;
		min = 0;
		for (int r = 0; r < rows; r++) {
			double p_state = (double) stateCount[states[r][agentNumber]] / (double) observations;
			// Entropy takes the negative log:
			localEntropy[r] = - Math.log(p_state) / log_2;
			average += localEntropy[r];
			if (localEntropy[r] > max) {
				max = localEntropy[r];
			} else if (localEntropy[r] < min) {
				min = localEntropy[r];
			}
		}
		average = average/(double) (rows);
		
		return localEntropy;
		
	}
	
	@Override
	public double[] computeLocalFromPreviousObservations(int states[][][], int agentIndex1, int agentIndex2){
		int timeSteps = states.length;
		//int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localEntropy = new double[timeSteps];
		average = 0;
		max = 0;
		min = 0;
		for (int r = 0; r < timeSteps; r++) {
			double p_state = (double) stateCount[states[r][agentIndex1][agentIndex2]] / (double) observations;
			// Entropy takes the negative log:
			localEntropy[r] = - Math.log(p_state) / log_2;
			average += localEntropy[r];
			if (localEntropy[r] > max) {
				max = localEntropy[r];
			} else if (localEntropy[r] < min) {
				min = localEntropy[r];
			}
		}
		average = average/(double) (timeSteps);
		
		return localEntropy;
		
	}

	@Override
	public final double[] computeLocal(int states[]) {
		
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}

	@Override
	public final double[][] computeLocal(int states[][]) {
		
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}

	@Override
	public final double[][][] computeLocal(int states[][][]) {
		
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}

	@Override
	public final double computeAverageLocal(int states[]) {
		
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}

	@Override
	public final double computeAverageLocal(int states[][]) {
		
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}

	@Override
	public final double computeAverageLocal(int states[][][]) {
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}

	@Override
	public final double[] computeLocal(int states[][], int col) {
		initialise();
		addObservations(states, col);
		return computeLocalFromPreviousObservations(states, col);
	}

	@Override
	public final double[] computeLocal(int states[][][],
			int agentIndex1, int agentIndex2) {
		initialise();
		addObservations(states, agentIndex1, agentIndex2);
		return computeLocalFromPreviousObservations(states, agentIndex1, agentIndex2);
	}

	@Override
	public final double computeAverageLocal(int states[][], int col) {
		initialise();
		addObservations(states, col);
		return computeAverageLocalOfObservations();
	}

	@Override
	public final double computeAverageLocal(int states[][][], int agentIndex1, int agentIndex2) {
		initialise();
		addObservations(states, agentIndex1, agentIndex2);
		return computeAverageLocalOfObservations();
	}
}
