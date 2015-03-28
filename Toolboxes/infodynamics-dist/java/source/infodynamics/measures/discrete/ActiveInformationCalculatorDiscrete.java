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
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

/**
 * <p>Active Information Storage  calculator for univariate discrete (int[]) data.
 * See definition of Active Information Storage (AIS) by Lizier et al. below.
 * Basically, AIS is the mutual information between the past <i>state</i>
 * of a time-series process <i>X</i> and its next value. The past <i>state</i> at time <code>n</code>
 * is represented by an embedding vector of <code>k</code> values from <code>X_n</code> backwards,
 * each separated by <code>\tau</code> steps, giving
 * <code><b>X^k_n</b> = [ X_{n-(k-1)\tau}, ... , X_{n-\tau}, X_n]</code>.
 * We call <code>k</code> the embedding dimension, and <code>\tau</code>
 * the embedding delay (only delay = 1 is implemented at the moment).
 * AIS is then the mutual information between <b>X^k_n</b> and X_{n+1}.</p>
 *
 * <p>Usage of the class is intended to follow this paradigm:</p>
 * <ol>
 * 		<li>Construct the calculator: {@link #ActiveInformationCalculatorDiscrete(int, int)};</li>
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
 * 	<li>J.T. Lizier, M. Prokopenko and A.Y. Zomaya,
 * 		<a href="http://dx.doi.org/10.1016/j.ins.2012.04.016">
 * 		"Local measures of information storage in complex distributed computation"</a>,
 * 		Information Sciences, vol. 208, pp. 39-54, 2012.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class ActiveInformationCalculatorDiscrete extends SingleAgentMeasureDiscreteInContextOfPastCalculator {

	/**
	 * User was formerly forced to create new instances through this factory method.
	 * Retained for backwards compatibility.
	 * 
	 * @param base
	 * @param history
	 * @deprecated
	 * @return
	 */
	public static ActiveInformationCalculatorDiscrete newInstance(int base, int history) {
		return new ActiveInformationCalculatorDiscrete(base, history);
	}
	
	/**
	 * Construct a new instance
	 * 
	 * @param base number of symbols for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param history embedded history length of the destination to condition on -
	 *        this is k in Schreiber's notation.
	 */
	public ActiveInformationCalculatorDiscrete(int base, int history) {
		super(base, history);
	}

	@Override
	public void addObservations(int states[]) {
		int timeSteps = states.length;
		// increment the count of observations:
		observations += (timeSteps - k); 
		
		// Initialise and store the current previous value for each column
		int prevVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p];
		}
		
		// 1. Count the tuples observed
		int nextVal;
		for (int t = k; t < timeSteps; t++) {
			// Add to the count for this particular transition:
			nextVal = states[t];
			nextPastCount[nextVal][prevVal]++;
			pastCount[prevVal]++;
			nextCount[nextVal]++;					
			// Update the previous value:
			prevVal -= maxShiftedValue[states[t-k]];
			prevVal *= base;
			prevVal += states[t];
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
				nextCount[nextVal]++;					
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
					nextCount[nextVal]++;					
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
			nextCount[nextVal]++;					
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
		
		// Initialise and store the current previous value for each column
		int prevVal = 0; 
		prevVal = 0;
		for (int p = 0; p < k; p++) {
			prevVal *= base;
			prevVal += states[p][agentIndex1][agentIndex2];
		}
		
		// 1. Count the tuples observed
		int nextVal;
		for (int t = k; t < timeSteps; t++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			nextVal = states[t][agentIndex1][agentIndex2];
			nextPastCount[nextVal][prevVal]++;
			pastCount[prevVal]++;
			nextCount[nextVal]++;					
			// Update the previous value:
			prevVal -= maxShiftedValue[states[t-k][agentIndex1][agentIndex2]];
			prevVal *= base;
			prevVal += states[t][agentIndex1][agentIndex2];
		}
	}

	@Override
	public double computeAverageLocalOfObservations() {
		double mi = 0.0;
		double miCont = 0.0;

		max = 0;
		min = 0;
		for (int nextVal = 0; nextVal < base; nextVal++) {
			// compute p_next
			double p_next = (double) nextCount[nextVal] / (double) observations;
			for (int prevVal = 0; prevVal < base_power_k; prevVal++) {
				// compute p_prev
				double p_prev = (double) pastCount[prevVal] / (double) observations;
				// compute p(prev, next)
				double p_joint = (double) nextPastCount[nextVal][prevVal] / (double) observations;
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
	 * Returns the average local entropy rate from
	 *  the observed values which have been passed in previously. 
	 *  
	 * @return average entropy rate
	 */
	public double computeAverageLocalEntropyRateOfObservations() {
		double entRate = 0.0;
		double entRateCont = 0.0;

		for (int nextVal = 0; nextVal < base; nextVal++) {
			for (int prevVal = 0; prevVal < base_power_k; prevVal++) {
				// compute p_prev
				double p_prev = (double) pastCount[prevVal] / (double) observations;
				// compute p(prev, next)
				double p_joint = (double) nextPastCount[nextVal][prevVal] / (double) observations;
				// Compute entropy rate contribution:
				if (p_joint > 0.0) {
					double logTerm = p_joint / p_prev;
					// Entropy rate takes the negative log:
					double localValue = - Math.log(logTerm) / log_base;
					entRateCont = p_joint * localValue;
				} else {
					entRateCont = 0.0;
				}
				entRate += entRateCont;
			}
		}
		
		return entRate;
	}
	
	/**
	 * Computes local active info storage for the given (single)
	 *  specific values.
	 *  
	 * See {@link TransferEntropyCalculatorDiscrete#getPastCount(int)} for how the
	 *  joint embedded values representing the past are calculated.
	 * 
	 * @param next next value of the variable
	 * @param past int representing the joint state of the past of the variable x[n]^k
	 * @return local active info storage value
	 */
	public double computeLocalFromPreviousObservations(int next, int past){
		double logTerm = ( (double) nextPastCount[next][past] ) /
		  ( (double) nextCount[next] *
			(double) pastCount[past] );
		logTerm *= (double) observations;
		return Math.log(logTerm) / log_base;
	}

	@Override
	public double[] computeLocalFromPreviousObservations(int states[]){
		int timeSteps = states.length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localActive = new double[timeSteps];
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
		for (int t = k; t < timeSteps; t++) {
			nextVal = states[t];
			logTerm = ( (double) nextPastCount[nextVal][prevVal] ) /
			  		  ( (double) nextCount[nextVal] *
			  			(double) pastCount[prevVal] );
			// Now account for the fact that we've
			//  just used counts rather than probabilities,
			//  and we've got two counts on the bottom
			//  but one count on the top:
			logTerm *= (double) observations;
			localActive[t] = Math.log(logTerm) / log_base;
			average += localActive[t];
			if (localActive[t] > max) {
				max = localActive[t];
			} else if (localActive[t] < min) {
				min = localActive[t];
			}
			// Update the previous value:
			prevVal -= maxShiftedValue[states[t-k]];
			prevVal *= base;
			prevVal += states[t];
		}
		average = average/(double) (timeSteps - k);
		
		return localActive;
		
	}

	@Override
	public double[][] computeLocalFromPreviousObservations(int states[][]){
		int rows = states.length;
		int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][] localActive = new double[rows][columns];
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
				  		  ( (double) nextCount[nextVal] *
				  			(double) pastCount[prevVal[c]] );
				// Now account for the fact that we've
				//  just used counts rather than probabilities,
				//  and we've got two counts on the bottom
				//  but one count on the top:
				logTerm *= (double) observations;
				localActive[r][c] = Math.log(logTerm) / log_base;
				average += localActive[r][c];
				if (localActive[r][c] > max) {
					max = localActive[r][c];
				} else if (localActive[r][c] < min) {
					min = localActive[r][c];
				}
				// Update the previous value:
				prevVal[c] -= maxShiftedValue[states[r-k][c]];
				prevVal[c] *= base;
				prevVal[c] += states[r][c];
			}
		}
		average = average/(double) (columns * (rows - k));
		
		return localActive;
		
	}
	
	@Override
	public double[][][] computeLocalFromPreviousObservations(int states[][][]){
		int timeSteps = states.length;
		int agentRows = states[0].length;
		int agentColumns = states[0][0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][][] localActive = new double[timeSteps][agentRows][agentColumns];
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
					  		  ( (double) nextCount[nextVal] *
					  			(double) pastCount[prevVal[r][c]] );
					// Now account for the fact that we've
					//  just used counts rather than probabilities,
					//  and we've got two counts on the bottom
					//  but one count on the top:
					logTerm *= (double) observations;
					localActive[t][r][c] = Math.log(logTerm) / log_base;
					average += localActive[t][r][c];
					if (localActive[t][r][c] > max) {
						max = localActive[t][r][c];
					} else if (localActive[t][r][c] < min) {
						min = localActive[t][r][c];
					}
					// Update the previous value:
					prevVal[r][c] -= maxShiftedValue[states[t-k][r][c]];
					prevVal[r][c] *= base;
					prevVal[r][c] += states[t][r][c];
				}
			}
		}
		average = average/(double) (agentRows * agentColumns * (timeSteps - k));
		
		return localActive;
	}

	@Override
	public double[] computeLocalFromPreviousObservations(int states[][], int col){
		int rows = states.length;
		//int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localActive = new double[rows];
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
			  		  ( (double) nextCount[nextVal] *
			  			(double) pastCount[prevVal] );
			// Now account for the fact that we've
			//  just used counts rather than probabilities,
			//  and we've got two counts on the bottom
			//  but one count on the top:
			logTerm *= (double) observations;
			localActive[r] = Math.log(logTerm) / log_base;
			average += localActive[r];
			if (localActive[r] > max) {
				max = localActive[r];
			} else if (localActive[r] < min) {
				min = localActive[r];
			}
			// Update the previous value:
			prevVal -= maxShiftedValue[states[r-k][col]];
			prevVal *= base;
			prevVal += states[r][col];
		}
		average = average/(double) (rows - k);
		
		return localActive;
		
	}
	
	@Override
	public double[] computeLocalFromPreviousObservations(int states[][][],
			int agentIndex1, int agentIndex2){
		int timeSteps = states.length;
		//int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localActive = new double[timeSteps];
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
			  		  ( (double) nextCount[nextVal] *
			  			(double) pastCount[prevVal] );
			// Now account for the fact that we've
			//  just used counts rather than probabilities,
			//  and we've got two counts on the bottom
			//  but one count on the top:
			logTerm *= (double) observations;
			localActive[t] = Math.log(logTerm) / log_base;
			average += localActive[t];
			if (localActive[t] > max) {
				max = localActive[t];
			} else if (localActive[t] < min) {
				min = localActive[t];
			}
			// Update the previous value:
			prevVal -= maxShiftedValue[states[t-k][agentIndex1][agentIndex2]];
			prevVal *= base;
			prevVal += states[t][agentIndex1][agentIndex2];
		}
		average = average/(double) (timeSteps - k);
		
		return localActive;
		
	}

	/**
	 * Generate a bootstrapped distribution of what the AIS would look like,
	 * under a null hypothesis that the past <code>k</code> values of our
	 * samples had no relation to the next value.
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for 
	 * an MI (like the AIS).
	 * </p>
	 * 
	 * <p>Note that if several disjoint time-series have been added 
	 * as observations using {@link #addObservations(int[], int[])} etc.,
	 * then these separate "trials" will be mixed up in the generation
	 * of surrogates here.</p>
	 * 
	 * <p>This method (in contrast to {@link #computeSignificance(int[][])})
	 * creates <i>random</i> shufflings of the next values for the surrogate AIS
	 * calculations.</p>
	 * 
	 * @param numPermutationsToCheck number of surrogate samples to bootstrap
	 *  to generate the distribution.
	 * @return the distribution of AIS scores under this null hypothesis.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) {
		RandomGenerator rg = new RandomGenerator();
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(observations, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}
	
	/**
	 * Generate a bootstrapped distribution of what the AIS would look like,
	 * under a null hypothesis that the previous <code>k</code> values of our
	 * samples had no relation to the next value in the time-series.
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for AIS 
	 * as a mutual information. Basically, the marginal PDFs
	 * of the past <code>k</code> values, and that of the next value, 
	 * are preserved, while their joint PDF is destroyed, and the 
	 * distribution of AIS under these conditions is generated.</p>
	 * 
	 * <p>Note that if several disjoint time-series have been added 
	 * as observations using {@link #addObservations(double[])} etc.,
	 * then these separate "trials" will be mixed up in the generation
	 * of surrogates here.</p>
	 * 
	 * <p>This method (in contrast to {@link #computeSignificance(int)})
	 * allows the user to specify how to construct the surrogates,
	 * such that repeatable results may be obtained.</p>
	 * 
	 * @param newOrderings a specification of how to shuffle the next values
	 *  to create the surrogates to generate the distribution with. The first
	 *  index is the permutation number (i.e. newOrderings.length is the number
	 *  of surrogate samples we use to bootstrap to generate the distribution here.)
	 *  Each array newOrderings[i] should be an array of length N (where
	 *  would be the value returned by {@link #getNumObservations()}),
	 *  containing a permutation of the values in 0..(N-1).
	 * @return the distribution of AIS scores under this null hypothesis.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception where the length of each permutation in newOrderings
	 *   is not equal to the number N samples that were previously supplied.
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) {
		double actualMI = computeAverageLocalOfObservations();
		
		int numPermutationsToCheck = newOrderings.length;
		
		// Reconstruct the values of the previous and next variables (not necessarily in order)
		int[] prevValues = new int[observations];
		int[] nextValues = new int[observations];
		int t_prev = 0;
		int t_next = 0;
		for (int prevVal = 0; prevVal < pastCount.length; prevVal++) {
			int numberOfSamplesPrev = pastCount[prevVal];
			MatrixUtils.fill(prevValues, prevVal, t_prev, numberOfSamplesPrev);
			t_prev += numberOfSamplesPrev;
		}
		for (int nextVal = 0; nextVal < base; nextVal++) {
			int numberOfSamplesNext = nextCount[nextVal];
			MatrixUtils.fill(nextValues, nextVal, t_next, numberOfSamplesNext);
			t_next += numberOfSamplesNext;
		}
		
		ActiveInformationCalculatorDiscrete ais2;
		ais2 = new ActiveInformationCalculatorDiscrete(base, k);
		ais2.initialise();
		ais2.observations = observations;
		ais2.pastCount = pastCount;
		ais2.nextCount = nextCount;
		int countWhereMIIsMoreSignificantThanOriginal = 0;
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		for (int p = 0; p < numPermutationsToCheck; p++) {
			// Generate a new re-ordered data set for the next variable
			int[] newDataNext = MatrixUtils.extractSelectedTimePoints(nextValues, newOrderings[p]);
			// compute the joint probability distribution
			MatrixUtils.fill(ais2.nextPastCount, 0);
			for (int t = 0; t < observations; t++) {
				ais2.nextPastCount[newDataNext[t]][prevValues[t]]++;
			}
			// And get an MI value for this realisation:
			double newMI = ais2.computeAverageLocalOfObservations();
			measDistribution.distribution[p] = newMI;
			if (newMI >= actualMI) {
				countWhereMIIsMoreSignificantThanOriginal++;
			}

		}
		
		// And return the significance
		measDistribution.pValue = (double) countWhereMIIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualMI;
		return measDistribution;
	}

	/**
	 * Debug method to write the current probability distribution functions 
	 *  
	 * @return
	 */
	public void writePdfs() {
		double mi = 0.0;
		double miCont = 0.0;

		System.out.println("nextVal p(next) prevVal p(prev) p(joint) logTerm localVal");
		for (int nextVal = 0; nextVal < base; nextVal++) {
			// compute p_next
			double p_next = (double) nextCount[nextVal] / (double) observations;
			for (int prevVal = 0; prevVal < base_power_k; prevVal++) {
				// compute p_prev
				double p_prev = (double) pastCount[prevVal] / (double) observations;
				// compute p(prev, next)
				double p_joint = (double) nextPastCount[nextVal][prevVal] / (double) observations;
				// Compute MI contribution:
				if (p_joint * p_next * p_prev > 0.0) {
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
	 * Utility function to compute a unique number to represent the
	 * combined past values of x up to and including time step t:
	 *  (i.e. (x_{t-k+1}, ... ,x_{t-1},x_{t}))
	 * 
	 * See {@link TransferEntropyCalculatorDiscrete#getPastCount(int)} for
	 *  how the joint value representing the past is calculated.
	 * 
	 * @param x time series
	 * @param t time step at which to compute the combined past
	 * @return an int representing the joint state of the past of x, x[t]^k
	 * 
	 */
	public int computePastValue(int[] x, int t) {
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += x[t - k + 1 + p];
		}
		return pastVal;
	}

	/**
	 * Utility function to compute a unique number to represent the
	 * combined past values of x (which is a column in data)
	 * up to and including time step t:
	 *  (i.e. (x_{t-k+1}, ... ,x_{t-1},x_{t}))
	 * 
	 * See {@link TransferEntropyCalculatorDiscrete#getPastCount(int)} for
	 *  how the joint value representing the past is calculated.
	 *  
	 * @param data 2D time series, first index is time,
	 *    second is variable number
	 * @param column column of data which is variable x
	 * @param t time step at which to compute the combined past
	 * @return an int representing the joint state of the past of x, x[t]^k
	 */
	public int computePastValue(int[][] data, int column, int t) {
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += data[t - k + 1 + p][column];
		}
		return pastVal;
	}

	/**
	 * Utility function to compute a unique number to represent the
	 * combined past values of x (which is a variable in data)
	 * up to and including time step t:
	 *  (i.e. (x_{t-k+1}, ... ,x_{t-1},x_{t}))
	 * 
	 * See {@link TransferEntropyCalculatorDiscrete#getPastCount(int)} for
	 *  how the joint value representing the past is calculated.
	 *  
	 * @param data 3D time series, first index is time,
	 *    second is variable row number, third is variable column number
	 * @param agentRow row of data for variable x
	 * @param agentColumn column of data for variable x
	 * @param t time step at which to compute the combined past
	 * @return an int representing the joint state of the past of x, x[t]^k
	 * Utility function to compute the combined past values of x up to and including time step t
	 *  (i.e. (x_{t-k+1}, ... ,x_{t-1},x_{t}))
	 */
	public int computePastValue(int[][][] data, int agentRow,
			int agentColumn, int t) {
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += data[t - k + 1 + p][agentRow][agentColumn];
		}
		return pastVal;
	}
}
