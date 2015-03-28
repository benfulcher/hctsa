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
 * A combined calculator for the active information,
 * entropy rate and entropy.
 * 
 * <p>This class is preliminary, so the Javadocs are incomplete --
 * please see {@link ActiveInformationCalculatorDiscrete},
 * {@link EntropyRateCalculatorDiscrete} and {@link EntropyCalculatorDiscrete}
 * for documentation on the corresponding functions
 * and typical usage pattern.
 * </p>
 * 
 * TODO Make this inherit from {@link SingleAgentMeasureDiscreteInContextOfPastCalculator}
 * like {@link ActiveInformationCalculatorDiscrete} and fix the Javadocs
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class CombinedActiveEntRateCalculatorDiscrete {

	private double averageActive = 0.0;
	private double maxActive = 0.0;
	private double minActive = 0.0;
	private double averageEntRate = 0.0;
	private double maxEntRate = 0.0;
	private double minEntRate = 0.0;
	private double averageEntropy = 0.0;
	private double maxEntropy = 0.0;
	private double minEntropy = 0.0;
	
	private int observations = 0;
	private int k = 0; // history length k. Need initialised to 0 for changedSizes
	private int base = 0; // number of individual states. Need initialised to 0 for changedSizes
	private int[][]	jointCount = null; // Count for (i[t+1], i[t]) tuples
	private int[] prevCount = null; // Count for i[t]		
	private int[] nextCount = null; // Count for i[t+1]

	// Space-time results (ST)
	// - for a homogeneous multi-agent system
	public class CombinedActiveEntRateLocalSTResults {
		public double[][] localActiveInfo;
		public double[][] localEntropyRate;
		public double[][] localEntropy;
	}
	
	public class CombinedActiveEntRateLocalResults {
		public double[] localActiveInfo;
		public double[] localEntropyRate;
		public double[] localEntropy;
	}

	public CombinedActiveEntRateCalculatorDiscrete() {
		super();
	}

	/**
	 * Initialise calculator, preparing to take observation sets in
	 * Should be called prior to any of the addObservations() methods.
	 * You can reinitialise without needing to create a new object.
	 *
	 */
	public void initialise(int history, int base){
		averageActive = 0.0;
		maxActive = 0.0;
		minActive = 0.0;
		observations = 0;
		
		boolean changedSizes = true;
		if ((this.base == base) && (this.k == history)) {
			changedSizes = false;
		}
		this.base = base;
		k = history;
		
		if (history < 1) {
			throw new RuntimeException("History k " + history + " is not >= 1 for ActiveInfo Calculator");
		}
		
		if (changedSizes) {
			// Create storage for counts of observations
			jointCount = new int[base][MathsUtils.power(base, history)];
			prevCount = new int[MathsUtils.power(base, history)];
			nextCount = new int[base];
		} else {
			// Just set counts to zeros without recreating the space
			MatrixUtils.fill(jointCount, 0);
			MatrixUtils.fill(prevCount, 0);
			MatrixUtils.fill(nextCount, 0);
		}
	}
		
	
	/**
 	 * Add observations in to our estimates of the pdfs.
 	 * This call suitable only for homogeneous agents, as all
 	 *  agents will contribute to single pdfs.
	 *
	 * @param states
	 */
	public void addObservations(int states[][]) {
		int rows = states.length;
		int columns = states[0].length;
		// increment the count of observations:
		observations += (rows - k)*columns; 
		
		// 1. Count the tuples observed
		int prevVal, nextVal;
		for (int r = k; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				// Add to the count for this particular transition:
				// (cell's assigned as above)
				nextVal = states[r][c];
				prevVal = 0;
				int multiplier = 1;
				for (int p = 1; p <= k; p++) {
					prevVal += states[r-p][c] * multiplier;			
					multiplier *= base;
				}
				jointCount[nextVal][prevVal]++;
				prevCount[prevVal]++;
				nextCount[nextVal]++;					
			}
		}		
	}
	
	/**
 	 * Add observations for a single agent of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to addObservations(int states[][])
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states
	 */
	public void addObservations(int states[][], int col) {
		int rows = states.length;
		// increment the count of observations:
		observations += (rows - k); 
		
		// 1. Count the tuples observed
		int prevVal, nextVal;
		for (int r = k; r < rows; r++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			nextVal = states[r][col];
			prevVal = 0;
			int multiplier = 1;
			for (int p = 1; p <= k; p++) {
				prevVal += states[r-p][col] * multiplier;			
				multiplier *= base;
			}
			jointCount[nextVal][prevVal]++;
			prevCount[prevVal]++;
			nextCount[nextVal]++;					
		}
	}

	/**
	 * Computes the average local values from
	 *  the observed values which have been passed in previously. 
	 * Access the averages, mins and maxes from the accessor methods.
	 *  
	 * @return
	 */
	public void computeAverageLocalOfObservations() {
		double activeCont, entropyCont;

		resetOverallStats();
		double localActiveValue, localEntropyValue, localEntRateValue, logTerm;
		
		for (int nextVal = 0; nextVal < base; nextVal++) {
			// compute p_next
			double p_next = (double) nextCount[nextVal] / (double) observations;
			// ** ENTROPY **
			if (p_next > 0.0) {
				// Entropy takes the negative log:
				localEntropyValue = - Math.log(p_next) / Math.log(base);
				entropyCont = p_next * localEntropyValue;
				if (localEntropyValue > maxEntropy) {
					maxEntropy = localEntropyValue;
				} else if (localEntropyValue < minEntropy) {
					minEntropy = localEntropyValue;
				}
			} else {
				localEntropyValue = 0.0;
				entropyCont = 0.0;
				continue; // no point computing ent rate and active info, will be zeros
			}
			averageEntropy += entropyCont;
			for (int prevVal = 0; prevVal < MathsUtils.power(base, k); prevVal++) {
				// compute p_prev
				double p_prev = (double) prevCount[prevVal] / (double) observations;
				// compute p(prev, next)
				double p_joint = (double) jointCount[nextVal][prevVal] / (double) observations;
				if (p_joint > 0.0) {
					// ** ACTIVE INFO **
					logTerm = p_joint / (p_next * p_prev);
					localActiveValue = Math.log(logTerm) / Math.log(base);
					activeCont = p_joint * localActiveValue;
					if (localActiveValue > maxActive) {
						maxActive = localActiveValue;
					} else if (localActiveValue < minActive) {
						minActive = localActiveValue;
					}
				} else {
					localActiveValue = 0.0;
					activeCont = 0.0;
				}
				averageActive += activeCont;
				// ** ENTROPY RATE ** = ENTROPY - ACTIVE
				localEntRateValue = localEntropyValue - localActiveValue;
				if (localEntRateValue > maxEntRate) {
					maxEntRate = localEntRateValue;
				} else if (localEntRateValue < maxEntRate) {
					maxEntRate = localEntRateValue;
				}
			}
		}
		averageEntRate = averageEntropy - averageActive;
		return;
	}
	
	/**
	 * Computes local values for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 * This method to be used for homogeneous agents only
	 *  
	 * @param states
	 * @return
	 */
	public CombinedActiveEntRateLocalSTResults computeLocalFromPreviousObservations(int states[][]){
		int rows = states.length;
		int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][] localActive = new double[rows][columns];
		double[][] localEntRate = new double[rows][columns];
		double[][] localEntropy = new double[rows][columns];
		
		resetOverallStats();
		
		int prevVal, nextVal;
		double logTerm;
		for (int r = k; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				nextVal = states[r][c];
				prevVal = 0;
				int multiplier = 1;
				for (int p = 1; p <= k; p++) {
					prevVal += states[r-p][c] * multiplier;			
					multiplier *= base;
				}
				// ** ENTROPY **
				double p_next = (double) nextCount[nextVal] / (double) observations;
				// Entropy takes the negative log:
				localEntropy[r][c] = - Math.log(p_next) / Math.log(base);
				averageEntropy += localEntropy[r][c];
				if (localEntropy[r][c] > maxEntropy) {
					maxEntropy = localEntropy[r][c];
				} else if (localEntropy[r][c] < minEntropy) {
					minEntropy = localEntropy[r][c];
				}
				// ** ACTIVE INFO **
				logTerm = ( (double) jointCount[nextVal][prevVal] ) /
				  		  ( (double) nextCount[nextVal] *
				  			(double) prevCount[prevVal] );
				// Now account for the fact that we've
				//  just used counts rather than probabilities,
				//  and we've got two counts on the bottom
				//  but one count on the top:
				logTerm *= (double) observations;
				localActive[r][c] = Math.log(logTerm) / Math.log(base);
				averageActive += localActive[r][c];
				if (localActive[r][c] > maxActive) {
					maxActive = localActive[r][c];
				} else if (localActive[r][c] < minActive) {
					minActive = localActive[r][c];
				}
				// ** ENTROPY RATE **
				localEntRate[r][c] = localEntropy[r][c] - localActive[r][c];
				if (localEntRate[r][c] > maxEntRate) {
					maxEntRate = localEntRate[r][c];
				} else if (localEntRate[r][c] < minEntRate) {
					minEntRate = localEntRate[r][c];
				}
			}
		}
		averageActive = averageActive/(double) (columns * (rows - k));
		averageEntropy = averageEntropy/(double) (columns * (rows - k));
		averageEntRate = averageEntropy - averageActive;
		
		// Package results ready for return
		CombinedActiveEntRateLocalSTResults results = new CombinedActiveEntRateLocalSTResults();
		results.localActiveInfo = localActive;
		results.localEntropyRate = localEntRate;
		results.localEntropy = localEntropy;
		
		return results;
	}
	
	/**
	 * Computes local values for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 * This method is suitable for heterogeneous agents
	 *  
	 * @param states
	 * @return
	 */
	public CombinedActiveEntRateLocalResults computeLocalFromPreviousObservations(int states[][], int col){
		int rows = states.length;
		//int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localActive = new double[rows];
		double[] localEntRate = new double[rows];
		double[] localEntropy = new double[rows];

		resetOverallStats();
		
		int prevVal, nextVal;
		double logTerm = 0.0;
		for (int r = k; r < rows; r++) {
			nextVal = states[r][col];
			prevVal = 0;
			int multiplier = 1;
			for (int p = 1; p <= k; p++) {
				prevVal += states[r-p][col] * multiplier;			
				multiplier *= base;
			}
			// ** ENTROPY **
			double p_next = (double) nextCount[nextVal] / (double) observations;
			// Entropy takes the negative log:
			localEntropy[r] = - Math.log(p_next) / Math.log(base);
			averageEntropy += localEntropy[r];
			if (localEntropy[r] > maxEntropy) {
				maxEntropy = localEntropy[r];
			} else if (localEntropy[r] < minEntropy) {
				minEntropy = localEntropy[r];
			}
			// ** ACTIVE INFO **
			logTerm = ( (double) jointCount[nextVal][prevVal] ) /
			  		  ( (double) nextCount[nextVal] *
			  			(double) prevCount[prevVal] );
			// Now account for the fact that we've
			//  just used counts rather than probabilities,
			//  and we've got two counts on the bottom
			//  but one count on the top:
			logTerm *= (double) observations;
			localActive[r] = Math.log(logTerm) / Math.log(base);
			averageActive += localActive[r];
			if (localActive[r] > maxActive) {
				maxActive = localActive[r];
			} else if (localActive[r] < minActive) {
				minActive = localActive[r];
			}
			// ** ENTROPY RATE **
			localEntRate[r] = localEntropy[r] - localActive[r];
			if (localEntRate[r] > maxEntRate) {
				maxEntRate = localEntRate[r];
			} else if (localEntRate[r] < minEntRate) {
				minEntRate = localEntRate[r];
			}
		}
		averageActive = averageActive/(double) (rows - k);
		averageEntropy = averageEntropy/(double) (rows - k);
		averageEntRate = averageEntropy - averageActive;

		// Package results ready for return
		CombinedActiveEntRateLocalResults results = new CombinedActiveEntRateLocalResults();
		results.localActiveInfo = localActive;
		results.localEntropyRate = localEntRate;
		results.localEntropy = localEntropy;

		return results;
		
	}
	
	/**
	 * Standalone routine to 
	 * compute local values across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * 
	 * @param history - parameter k
	 * @param base - base of the states
	 * @param states - 2D array of states
	 * @return
	 */
	public CombinedActiveEntRateLocalSTResults computeLocal(int history, int base, int states[][]) {
		
		initialise(history, base);
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}
	
	/**
	 * Standalone routine to 
	 * compute average local values across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return the average
	 * This method to be called for homogeneous agents only
	 * 
	 * @param history - parameter k
	 * @param base - base of the states
	 * @param states - 2D array of states
	 * @return
	 */
	public void computeAverageLocal(int history, int base, int states[][]) {
		
		initialise(history, base);
		addObservations(states);
		computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute local values for one agent in a 2D spatiotemporal
	 *  array of the states of agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * This method should be used for heterogeneous agents
	 * 
	 * @param history - parameter k
	 * @param base - base of the states
	 * @param states - 2D array of states
	 * @param col - column number of the agent in the states array
	 * @return
	 */
	public CombinedActiveEntRateLocalResults computeLocal(int history, int base, int states[][], int col) {
		
		initialise(history, base);
		addObservations(states, col);
		return computeLocalFromPreviousObservations(states, col);
	}
	
	/**
	 * Standalone routine to 
	 * compute average local values
	 * for a single agent
	 * Returns the average
	 * This method suitable for heterogeneous agents
	 * 
	 * @param history - parameter k
	 * @param base - base of the states
	 * @param states - 2D array of states
	 * @param col - column number of the agent in the states array
	 * @return
	 */
	public void computeAverageLocal(int history, int base, int states[][], int col) {
		
		initialise(history, base);
		addObservations(states, col);
		computeAverageLocalOfObservations();
	}

	private void resetOverallStats() {
		averageActive = 0;
		maxActive = 0;
		minActive = 0;
		averageEntRate = 0;
		maxEntRate = 0;
		minEntRate = 0;
		averageEntropy = 0;
		maxEntropy = 0;
		minEntropy = 0;
	}
	
	/*
	 * Accessors for last active information computation
	 */
	public double getLastAverageActive() {
		return averageActive;
	}
	public double getLastMaxActive() {
		return maxActive;
	}
	public double getLastMinActive() {
		return minActive;
	}

	/*
	 * Accessors for last entropy rate computation
	 */
	public double getLastAverageEntRate() {
		return averageEntRate;
	}
	public double getLastMaxEntRate() {
		return maxEntRate;
	}
	public double getLastMinEntRate() {
		return minEntRate;
	}

	/*
	 * Accessors for last entropy computation
	 */
	public double getLastAverageEntropy() {
		return averageEntropy;
	}
	public double getLastMaxEntropy() {
		return maxEntropy;
	}
	public double getLastMinEntropy() {
		return minEntropy;
	}
}
