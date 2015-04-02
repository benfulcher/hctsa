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

import infodynamics.measures.continuous.MultiInfoCalculator;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

import java.util.Properties;

/**
 * Implements <b>separable information</b> (see Lizier et al, 2010, below).
 * 
 * Separable information is the sum of active information and apparent transfer entropy from every
 *  causal information contributor.
 *
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * <ol>
 * 		<li>Construct the calculator via
 * 			{@link #SeparableInfoCalculatorDiscrete(int, int, int)};</li>
 *		<li>Initialise the calculator using
 *			{@link #initialise()};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using one or more calls to
 * 			the set of {@link #addObservations(int[][], int[])} methods, then</li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average TE: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local TE values for these samples: {@link #computeLocalOfPreviousObservations()}</li>
 * 				<li>local TE values for a specific set of samples: e.g.
 * 				{@link #computeLocalFromPreviousObservations(int[][], int[])} etc.</li>
 * 			</ul>
 * 		</li>
 * 		<li>As an alternative to steps 3 and 4, the user may undertake
 * 			standalone computation from a single set of observations, via
 *  		e.g.: {@link #computeLocal(int[][], int[])} or
 *  		{@link #computeAverageLocal(int[][], int[])}.</li>
 * 		<li>
 * 		Return to step 2 to re-use the calculator on a new data set.
 * 		</li>
 * 	</ol>
 * 
 * <p>The causal information contributors (either their offsets or their absolute column numbers)
 *  should be supplied in the same order in every method call, otherwise the answer supplied will
 *  be incorrect.</p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href=http://dx.doi.org/10.1063/1.3486801">
 *  "Information modification and particle collisions in distributed computation"</a>
 *  Chaos 20, 3, 037109 (2010).</li>
 * </ul>
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class SeparableInfoCalculatorDiscrete extends ContextOfPastMeasureCalculatorDiscrete {

	protected int numSources = 0;
	protected int base_power_sources = 0;
	protected int[][][][] sourceNumValueNextPastCount = null;	// count for (i-j[n],i[n+1],i[n]^k) tuples for each source
	protected int[][][] sourcesNextPastCount = null;			// count for ({i-j[n]},i[n+1],i[n]^k) tuples for set of all sources
	protected int[][][] sourceNumValuePastCount = null;			// count for (i-j[n],i[n]^k) tuples for each source
	// Extra items to measure:
	protected double avPositiveLocal = 0.0;
	protected double avNegativeLocal = 0.0;
	// Variables used during the computation of the averages that need storage due to recursion:
	protected double meanSqLocals;
	protected boolean periodicBoundaryConditions = true;
	// Track whether the user will want to ask for the coherence of computation or not.
	//  Use this because it takes extra time in the separable info calculation that would be
	//  wasted if there was no intention of asking for the coherence result.
	protected boolean computeMultiInfoCoherence = false;
	protected MultiInfoCalculator miCalc = null;
	
	// Indices for row and column agent indices in lists of source agents
	public static final int ROW_INDEX = 0;
	public static final int COLUMN_INDEX = 0;
	
	// Cut off for running separable info calculator by 
	//  adding active info and transfer entropies rather than
	//  directly multiplying logarithms. (Adding is slower but uses
	//  less memory.
	// Leaving this as public so that calling applications can 
	//  reset it if needs be.
	public static int MAX_CONFIGS_FOR_DIRECT_CALC = 1000000;
	// Or we could just force the use of the direct calculator
	public static boolean FORCE_DIRECT_CALC = false;
	
	/**
	 * User to create new instances through this factory method.
	 * This allows us to return an efficient calculator for
	 * particular situations without the user needing to have 
	 * knowledge of this.
	 * @param base
	 * @param history
	 * @param numInfoContributors
	 * 
	 * @return
	 */
	public static SeparableInfoCalculatorDiscrete
		newInstance(int base, int history, int numInfoContributors) {
		
		if (history < 1) {
			//TODO make this class compatible with k==0
			//  (low priority, not truly necessary) 
			throw new RuntimeException("This class does not currently " +
					"function with k < 1 (see ConditionalTransferEntropyCalculator " +
					"for how to implement this)");
		}
		
		if (!FORCE_DIRECT_CALC &&
				(MathsUtils.power(base,numInfoContributors + history + 1)
				> MAX_CONFIGS_FOR_DIRECT_CALC)) {
			return new SeparableInfoCalculatorDiscreteByAddition
				(base, history, numInfoContributors);
		} else {
			return new SeparableInfoCalculatorDiscrete
					(base, history, numInfoContributors);
		}
	}

	/**
	 * Construct an instance
	 * 
	 * @param base number of symbols for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param history embedded history length of the destination to condition on -
	 *        this is k in Schreiber's notation.
	 * @param numOtherInfoContributors number of information contributors
	 *   (other than the past of the destination
	 *   or the source) to condition on.
	 */
	public SeparableInfoCalculatorDiscrete
		(int base, int history, int numInfoContributors) {
		this(base, history, numInfoContributors, false);
	}

	/**
	 * Private method to implement the public constructor
	 * 
	 * @param base number of symbols for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param history embedded history length of the destination to condition on -
	 *        this is k in Schreiber's notation.
	 * @param numOtherInfoContributors number of information contributors
	 *   (other than the past of the destination
	 *   or the source) to condition on.
	 * @param dontCreateObsStorage indicates that storage for
	 *  observations should not be created.
	 */
	protected SeparableInfoCalculatorDiscrete
		(int base, int history, int numInfoContributors, boolean dontCreateObsStorage) {

		super(base, history, dontCreateObsStorage);
		numSources = numInfoContributors;
		base_power_sources = MathsUtils.power(base, numSources);
		if (numInfoContributors < 1) {
			throw new RuntimeException("Number of info contributors < 1 for SeparableInfoCalculator");
		}
		if (!dontCreateObsStorage) {
			// Create storage for extra counts of observations
			sourceNumValueNextPastCount = new int[numInfoContributors][base][base][base_power_k];
			sourcesNextPastCount = new int[base_power_sources][base][base_power_k];
			sourceNumValuePastCount = new int[numInfoContributors][base][base_power_k];
		}
	}

	@Override
	public void initialise(){
		super.initialise();
		
		if (!noObservationStorage) {
			MatrixUtils.fill(sourceNumValueNextPastCount, 0);
			MatrixUtils.fill(sourcesNextPastCount, 0);
			MatrixUtils.fill(sourceNumValuePastCount, 0);
		}
		
		if (computeMultiInfoCoherence) {
			// We will be looking at multi-information between active info and apparent TE
			//  from all sources
			miCalc.initialise(numSources + 1);
		}
	}
	
	/**
 	 * Add observations in to our estimates of the pdfs.
 	 * This call suitable only for homogeneous agents, as all
 	 *  agents will contribute to single pdfs, and all are assumed
 	 *  to have other info contributors at same offsets.
 	 * 
	 * @param states multivariate time series, 1st index is time, 2nd index is agent number
	 * @param offsetOfDestFromSources offsets of the destination *from* causal information contributors.
	 * 		  (i.e. an offset of 1 means the destination is one index larger, or one to the right,
	 * 			than the source). 
	 *        sourcesOffsets is permitted to include 0, it will be ignored.
	 */
	public void addObservations(int states[][], int offsetOfDestFromSources[]) {
		addObservations(states, offsetOfDestFromSources, false);
	}
	private void addObservations(int states[][], int offsetOfDestFromSources[], boolean cleanedSources) {
		
		int[] cleanedSourcesOffsets;
		if (cleanedSources) {
			cleanedSourcesOffsets = offsetOfDestFromSources;
		} else {
			cleanedSourcesOffsets = cleanOffsetOfDestFromSources(offsetOfDestFromSources);
			// This call made redundant by cleanOffsetSources:
			// confirmEnoughOffsetSources(sourcesOffsets, j);
		}
		
		int timeSteps = states.length;
		int numAgents = states[0].length;
		
		int minAgentOffset = MatrixUtils.min(cleanedSourcesOffsets);
		int maxAgentOffset = MatrixUtils.max(cleanedSourcesOffsets);
		int nonPeriodicStartAgent = Math.max(0, maxAgentOffset);
		int nonPeriodicEndAgent = numAgents - 1 + Math.min(0, minAgentOffset);
		// increment the count of observations:
		if (periodicBoundaryConditions) {
			observations += (timeSteps - k) * numAgents;
		} else {
			observations += (timeSteps - k) * (nonPeriodicEndAgent - nonPeriodicStartAgent + 1);
		}
		
		// Initialise and store the current previous value for each column
		int[] pastVal = new int[numAgents]; 
		for (int c = 0; c < numAgents; c++) {
			pastVal[c] = 0;
			for (int p = 0; p < k; p++) {
				pastVal[c] *= base;
				pastVal[c] += states[p][c];
			}
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal, jointSourcesVal;
		for (int t = k; t < timeSteps; t++) {
			for (int c = periodicBoundaryConditions ? 0 : nonPeriodicStartAgent;
			 c < (periodicBoundaryConditions ? numAgents : nonPeriodicEndAgent + 1);
			 c++) {
				// Add to the count for this particular transition:
				// (cell's assigned as above)
				destVal = states[t][c];
				nextPastCount[destVal][pastVal[c]]++;
				pastCount[pastVal[c]]++;
				nextCount[destVal]++;
				jointSourcesVal = 0;
				for (int sIndex = 0; sIndex < numSources; sIndex++) {
					sourceVal = states[t-1][(c-cleanedSourcesOffsets[sIndex]+numAgents) % numAgents];
					sourceNumValueNextPastCount[sIndex][sourceVal][destVal][pastVal[c]]++;
					sourceNumValuePastCount[sIndex][sourceVal][pastVal[c]]++;
					jointSourcesVal *= base;
					jointSourcesVal += sourceVal;
				}
				sourcesNextPastCount[jointSourcesVal][destVal][pastVal[c]]++;
				// Update the previous value:
				pastVal[c] -= maxShiftedValue[states[t-k][c]];
				pastVal[c] *= base;
				pastVal[c] += states[t][c];
			}
		}		
	}
	
	/**
 	 * Add observations in to our estimates of the pdfs.
 	 * This call suitable only for homogeneous agents, as all
 	 *  agents will contribute to single pdfs, and all are assumed
 	 *  to have other info contributors at same offsets.
 	 * 
	 * @param states multivariate time series 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @param offsetOfDestFromSources 2D offsets of the destination *from* causal information contributors.
	 * 		  (i.e. an offset of 1 means the destination is one index larger, or one to the right,
	 * 			than the source). 
	 *        offsetOfDestFromSources is permitted to include 0, it will be ignored.
	 */
	public void addObservations(int states[][][], int offsetOfDestFromSources[][]) {
		addObservations(states, offsetOfDestFromSources, false);
	}
	private void addObservations(int states[][][], int offsetOfDestFromSources[][], boolean cleanedSources) {
		
		int[][] cleanedSourcesOffsets;
		if (cleanedSources) {
			cleanedSourcesOffsets = offsetOfDestFromSources;
		} else {
			cleanedSourcesOffsets = cleanOffsetOfDestFromSources(offsetOfDestFromSources);
			// This call made redundant by cleanOffsetSources:
			// confirmEnoughOffsetSources(sourcesOffsets);
		}
				
		int timeSteps = states.length;
		if (timeSteps == 0) {
			return;
		}
		int agentRows = states[0].length;
		if (agentRows == 0) {
			return;
		}
		int agentColumns = states[0][0].length;
		int minRowOffset = MatrixUtils.min(cleanedSourcesOffsets, ROW_INDEX);
		int maxRowOffset = MatrixUtils.max(cleanedSourcesOffsets, ROW_INDEX);
		int nonPeriodicStartRow = Math.max(0, maxRowOffset);
		int nonPeriodicEndRow = agentRows - 1 + Math.min(0, minRowOffset);
		int minColumnOffset = MatrixUtils.min(cleanedSourcesOffsets, COLUMN_INDEX);
		int maxColumnOffset = MatrixUtils.max(cleanedSourcesOffsets, COLUMN_INDEX);
		int nonPeriodicStartColumn = Math.max(0, maxColumnOffset);
		int nonPeriodicEndColumn = agentColumns - 1 + Math.min(0, minColumnOffset);
		// increment the count of observations:
		if (periodicBoundaryConditions) {
			observations += (timeSteps - k) * agentRows * agentColumns;
		} else {
			observations += (timeSteps - k) * (nonPeriodicEndRow - nonPeriodicStartRow + 1) *
					(nonPeriodicEndColumn - nonPeriodicStartColumn + 1);
		}
		
		// Initialise and store the current previous value for each column
		int[][] pastVal = new int[agentRows][agentColumns];
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
				pastVal[r][c] = 0;
				for (int p = 0; p < k; p++) {
					pastVal[r][c] *= base;
					pastVal[r][c] += states[p][r][c];
				}
			}
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal, jointSourcesVal;
		for (int t = k; t < timeSteps; t++) {
			for (int r = periodicBoundaryConditions ? 0 : nonPeriodicStartRow;
					 r < (periodicBoundaryConditions ? agentRows : nonPeriodicEndRow + 1);
					 r++) {
				for (int c = periodicBoundaryConditions ? 0 : nonPeriodicStartColumn;
						c < (periodicBoundaryConditions ? agentColumns : nonPeriodicEndColumn + 1);
						c++) {
					// Add to the count for this particular transition:
					// (cell's assigned as above)
					destVal = states[t][r][c];
					nextPastCount[destVal][pastVal[r][c]]++;
					pastCount[pastVal[r][c]]++;
					nextCount[destVal]++;
					jointSourcesVal = 0;
					for (int sIndex = 0; sIndex < numSources; sIndex++) {
						// Can safely do mod operations here - if not periodic boundary conditions
						//  and this source would have bounced over boundary, we would have skipped
						//  over this destination earlier
						sourceVal = states[t-1]
						        [(r-cleanedSourcesOffsets[sIndex][ROW_INDEX]+agentRows) % agentRows]
						        [(c-cleanedSourcesOffsets[sIndex][COLUMN_INDEX]+agentColumns) % agentColumns];
						sourceNumValueNextPastCount[sIndex][sourceVal][destVal][pastVal[r][c]]++;
						sourceNumValuePastCount[sIndex][sourceVal][pastVal[r][c]]++;
						jointSourcesVal *= base;
						jointSourcesVal += sourceVal;
					}
					sourcesNextPastCount[jointSourcesVal][destVal][pastVal[r][c]]++;
					// Update the previous value:
					pastVal[r][c] -= maxShiftedValue[states[t-k][r][c]];
					pastVal[r][c] *= base;
					pastVal[r][c] += states[t][r][c];
				}
			}
		}		
	}

	/**
 	 * Add observations for a single destination of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to addObservations(int states[][])
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states multivariate time series, 1st index is time, 2nd index is agent number
	 * @param destCol the destination index
	 * @param sourcesAbsolute array of the source indices
	 */
	public void addObservations(int states[][], int destCol, int[] sourcesAbsolute) {
		addObservations(states, destCol, sourcesAbsolute, false);
	}
	private void addObservations(int states[][], int destCol, int[] sourcesAbsolute, boolean cleanedSources) {

		int[] cleanedSourcesAbsolute;
		if (cleanedSources) {
			cleanedSourcesAbsolute = sourcesAbsolute;
		} else {
			cleanedSourcesAbsolute = cleanAbsoluteSources(sourcesAbsolute, destCol);
			// This call made redundant by cleanAbsoluteSources:
			// confirmEnoughAbsoluteSources(sourcesAbsolute, destCol);
		}
		
		int rows = states.length;
		// increment the count of observations:
		observations += (rows - k); 
		
		// Initialise and store the current previous value for each column
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += states[p][destCol];
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal, jointSourcesVal;
		for (int r = k; r < rows; r++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			destVal = states[r][destCol];
			nextPastCount[destVal][pastVal]++;
			pastCount[pastVal]++;			
			nextCount[destVal]++;
			jointSourcesVal = 0;
			for (int sIndex = 0; sIndex < numSources; sIndex++) {
				sourceVal = states[r-1][cleanedSourcesAbsolute[sIndex]];
				sourceNumValueNextPastCount[sIndex][sourceVal][destVal][pastVal]++;
				sourceNumValuePastCount[sIndex][sourceVal][pastVal]++;					
				jointSourcesVal *= base;
				jointSourcesVal += sourceVal;
			}
			sourcesNextPastCount[jointSourcesVal][destVal][pastVal]++;
			// Update the previous value:
			pastVal -= maxShiftedValue[states[r-k][destCol]];
			pastVal *= base;
			pastVal += states[r][destCol];
		}
	}

	/**
 	 * Add observations for a single destination pair of the 
 	 * multivariate time series
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to addObservations(int states[][])
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states multivariate time series: 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @param destAgentRow the destination index
	 * @param destAgentColumn the destination index
	 * @param sourcesAbsolute array of the source indices
	 */
	public void addObservations(int states[][][], int destAgentRow, int destAgentColumn, int[][] sourcesAbsolute) {
		addObservations(states, destAgentRow, destAgentColumn, sourcesAbsolute, false);
	}
	private void addObservations(int states[][][], int destAgentRow, int destAgentColumn, int[][] sourcesAbsolute, boolean cleanedSources) {

		int[][] cleanedSourcesAbsolute;
		if (cleanedSources) {
			cleanedSourcesAbsolute = sourcesAbsolute;
		} else {
			cleanedSourcesAbsolute = cleanAbsoluteSources(sourcesAbsolute, destAgentRow, destAgentColumn);
			// This call made redundant by cleanAbsoluteSources:
			// confirmEnoughAbsoluteSources(sourcesAbsolute, destCol);
		}
		
		int timeSteps = states.length;
		// increment the count of observations:
		observations += (timeSteps - k); 
		
		// Initialise and store the current previous value for each column
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += states[p][destAgentRow][destAgentColumn];
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal, jointSourcesVal;
		for (int t = k; t < timeSteps; t++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			destVal = states[t][destAgentRow][destAgentColumn];
			nextPastCount[destVal][pastVal]++;
			pastCount[pastVal]++;			
			nextCount[destVal]++;
			jointSourcesVal = 0;
			for (int sIndex = 0; sIndex < numSources; sIndex++) {
				sourceVal = states[t-1][cleanedSourcesAbsolute[sIndex][ROW_INDEX]]
				                        [cleanedSourcesAbsolute[sIndex][COLUMN_INDEX]];
				sourceNumValueNextPastCount[sIndex][sourceVal][destVal][pastVal]++;
				sourceNumValuePastCount[sIndex][sourceVal][pastVal]++;					
				jointSourcesVal *= base;
				jointSourcesVal += sourceVal;
			}
			sourcesNextPastCount[jointSourcesVal][destVal][pastVal]++;
			// Update the previous value:
			pastVal -= maxShiftedValue[states[t-k][destAgentRow][destAgentColumn]];
			pastVal *= base;
			pastVal += states[t][destAgentRow][destAgentColumn];
		}
	}

	@Override
	public synchronized double computeAverageLocalOfObservations() {
		max = 0;
		min = 0;
		meanSqLocals = 0;
		average = 0.0;
		avPositiveLocal = 0.0;
		avNegativeLocal = 0.0;
		if (computeMultiInfoCoherence) {
			miCalc.startAddObservations();
		}
		
		// Create space for the joint source values to run through:
		int[] sourceValues = new int[numSources];
		// Call our recursive function, asking it to compute the average over all source values
		computeAverageLocalOfObservations(sourceValues, 0);
		// Close off the individual observations
		try {
			if (computeMultiInfoCoherence) {
				miCalc.finaliseAddObservations();
			}
		} catch (Exception e) {
			// an exception would only be thrown if we changed the number of causal contributors here
			//  which simply will not happen. Just in case it does, we'll throw a runtime exception
			throw new RuntimeException("Number of causal contributors changed from intialisation to calculation!");
		}
		std = Math.sqrt(meanSqLocals - average * average);
		return average;
	}
	
	/**
	 * Private utility function.
	 * 
	 * <p>Updates the average, max, min and meanSq of locals for the separable information
	 *  for the given source values in sourceValues up to the index indexToModify over
	 *  all possible source values after the index indexToModify onwards. Uses recursion
	 *  on increasing indexToModify.
	 * Designed to be called with indexToModify == 0 to compute average over all states.
	 * 
	 * @param sourceValues
	 * @param indexToModify
	 * @return
	 */
	protected void computeAverageLocalOfObservations(int[] sourceValues, int indexToModify) {

		if (indexToModify < sourceValues.length) {
			// Assign values to our variables and make the recursive call
			for (int s = 0; s < base; s++) {
				sourceValues[indexToModify] = s;
				computeAverageLocalOfObservations(sourceValues, indexToModify + 1);
			}
			return;
		}
		
		// else there were no more source values to assign, so we carry out the computation with the
		//  given source values
		
		double logTerm, localValue, sepCont;
		int jointSourcesVal, sourceVal;
		// Compute the joint source value first:
		jointSourcesVal = 0;
		for (int sIndex = 0; sIndex < numSources; sIndex++) {
			jointSourcesVal *= base;
			jointSourcesVal += sourceValues[sIndex];
		}
		
		// At this point, we have a value for every source and have computed the joint source
		// value. Now, let's look at varying the destination's past and next state, and taking
		// contributions for each tuple.
		double[] localActAndTes = new double[numSources + 1];
		for (int pastVal = 0; pastVal < base_power_k; pastVal++) {
			for (int destVal = 0; destVal < base; destVal++) {
				if (sourcesNextPastCount[jointSourcesVal][destVal][pastVal] != 0) {
					// Add in the local active information storage first:
					logTerm = ( (double) nextPastCount[destVal][pastVal] ) /
			  		  ( (double) nextCount[destVal] *
			  			(double) pastCount[pastVal] );
					// Now account for the fact that we've just used counts rather than probabilities,
					//  and we've got two counts on the bottom but one count on the top:
					logTerm *= (double) observations;
					if (computeMultiInfoCoherence) {
						// Keep the local active info.
						// Adding natural logs, since we're going to normalize these anyway
						localActAndTes[0] = Math.log(logTerm);
					}
					// Then add in the local transfer entropy for each source:
					for (int sIndex = 0; sIndex < numSources; sIndex++) {
						sourceVal = sourceValues[sIndex];
						localActAndTes[sIndex+1] = ((double) sourceNumValueNextPastCount[sIndex][sourceVal][destVal][pastVal] / (double) sourceNumValuePastCount[sIndex][sourceVal][pastVal]) /
					 		((double) nextPastCount[destVal][pastVal] / (double) pastCount[pastVal]);
						logTerm *= localActAndTes[sIndex+1];
						if (computeMultiInfoCoherence) {
							// Keep this local TE
							localActAndTes[sIndex+1] = Math.log(localActAndTes[sIndex+1]);
						}
					}
					// Add in these local active and TE values for the coherence calculation
					if (computeMultiInfoCoherence) {
						// Need to add this observation in once for every time it occurs
						// (since we're looping over possible tuples here rather than observations)
						for (int i = 0; i < sourcesNextPastCount[jointSourcesVal][destVal][pastVal]; i++) {
							miCalc.addObservation(localActAndTes);
						}
					}
					localValue = Math.log(logTerm) / log_2;
					sepCont = (double) sourcesNextPastCount[jointSourcesVal][destVal][pastVal] /
								(double) observations * localValue;
					average += sepCont;
					if (sepCont >= 0.0) {
						avPositiveLocal += sepCont;
					} else {
						avNegativeLocal += sepCont;
					}
					if (localValue > max) {
						max = localValue;
					} else if (localValue < min) {
						min = localValue;
					}
					// Add this contribution to the mean 
					//  of the squared local values
					meanSqLocals += sepCont * localValue;
				}
			}
		}
	}
	
	/**
	 * Computes local separable information for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method to be used for homogeneous agents only
	 *  since it assumes the source offsets are the same
	 *  for all destinations (and all are added to the PDFs)
	 *  
	 * @param states multivariate time series 1st index is time, 2nd is agent index
	 * @param offsetOfDestFromSources offsets of the destination *from* causal information contributors.
	 * 		  (i.e. an offset of 1 means the destination is one index larger, or one to the right,
	 * 			than the source). 
	 *        sourcesOffsets is permitted to include 0, it will be ignored.
	 * @return multivariate time series of local separable information,
	 *  indexed as per states.
	 */
	public double[][] computeLocalFromPreviousObservations
		(int states[][], int offsetOfDestFromSources[]){
		
		return computeLocalFromPreviousObservations(states, offsetOfDestFromSources, false);
	}
	private double[][] computeLocalFromPreviousObservations
		(int states[][], int offsetOfDestFromSources[], boolean cleanedOffsets){
		
		int[] cleanedOffsetOfDestFromSources;
		if (cleanedOffsets) {
			cleanedOffsetOfDestFromSources = offsetOfDestFromSources;
		} else {
			cleanedOffsetOfDestFromSources = cleanOffsetOfDestFromSources(offsetOfDestFromSources);
			// This call made redundant by cleanOffsetSources:
			// confirmEnoughOffsetSources(othersOffsets);
		}

		int timeSteps = states.length;
		int numAgents = states[0].length;
		int minAgentOffset = MatrixUtils.min(cleanedOffsetOfDestFromSources);
		int maxAgentOffset = MatrixUtils.max(cleanedOffsetOfDestFromSources);
		int nonPeriodicStartAgent = Math.max(0, maxAgentOffset);
		int nonPeriodicEndAgent = numAgents - 1 + Math.min(0, minAgentOffset);

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][] localSep = new double[timeSteps][numAgents];
		average = 0.0;
		avPositiveLocal = 0.0;
		avNegativeLocal = 0.0;
		max = 0.0;
		min = 0.0;

		// Initialise and store the current previous value for each column
		int[] pastVal = new int[numAgents]; 
		for (int c = 0; c < numAgents; c++) {
			pastVal[c] = 0;
			for (int p = 0; p < k; p++) {
				pastVal[c] *= base;
				pastVal[c] += states[p][c];
			}
		}
		
		// Make a vector of the active and TE values for the coherence computation
		if (computeMultiInfoCoherence) {
			miCalc.startAddObservations();
		}
		double[] localActAndTes = new double[numSources + 1];
		
		int destVal, sourceVal;
		double logTerm;
		for (int t = k; t < timeSteps; t++) {
			for (int c = periodicBoundaryConditions ? 0 : nonPeriodicStartAgent;
			 c < (periodicBoundaryConditions ? numAgents : nonPeriodicEndAgent + 1);
			 c++) {
				destVal = states[t][c];
				// Add in the local active information storage:
				logTerm = ( (double) nextPastCount[destVal][pastVal[c]] ) /
		  		  ( (double) nextCount[destVal] *
		  			(double) pastCount[pastVal[c]] );
				// Now account for the fact that we've just used counts rather than probabilities,
				//  and we've got two counts on the bottom but one count on the top:
				logTerm *= (double) observations;
				if (computeMultiInfoCoherence) {
					// Keep the local active info.
					// Adding natural logs, since we're going to normalize these anyway
					localActAndTes[0] = Math.log(logTerm);
				}
				for (int sIndex = 0; sIndex < numSources; sIndex++) {
					sourceVal = states[t-1][(c-cleanedOffsetOfDestFromSources[sIndex]+numAgents) % numAgents];
					localActAndTes[sIndex+1] = ((double) sourceNumValueNextPastCount[sIndex][sourceVal][destVal][pastVal[c]] / (double) sourceNumValuePastCount[sIndex][sourceVal][pastVal[c]]) /
			 			((double) nextPastCount[destVal][pastVal[c]] / (double) pastCount[pastVal[c]]);
					// Add in the local transfer entropy for each source:
					logTerm *= localActAndTes[sIndex+1];
					if (computeMultiInfoCoherence) {
						// Keep this local TE
						localActAndTes[sIndex+1] = Math.log(localActAndTes[sIndex+1]);
					}
				}
				// Add in these local active and TE values for the coherence calculation
				if (computeMultiInfoCoherence) {
					miCalc.addObservation(localActAndTes);
				}
				localSep[t][c] = Math.log(logTerm) / log_2;
				average += localSep[t][c];
				if (localSep[t][c] > 0.0) {
					avPositiveLocal += localSep[t][c];
				} else {
					avNegativeLocal += localSep[t][c];
				}
				if (localSep[t][c] > max) {
					max = localSep[t][c];
				} else if (localSep[t][c] < min) {
					min = localSep[t][c];
				}
				// Update the previous value:
				pastVal[c] -= maxShiftedValue[states[t-k][c]];
				pastVal[c] *= base;
				pastVal[c] += states[t][c];
			}
		}		
		if (periodicBoundaryConditions) {
			average /= (double) (numAgents * (timeSteps - k));
			avPositiveLocal /= (double) (numAgents * (timeSteps - k));
			avNegativeLocal /= (double) (numAgents * (timeSteps - k));
		} else {
			average /= (double) ((timeSteps - k) * (nonPeriodicEndAgent - nonPeriodicStartAgent + 1));
			avPositiveLocal /= (double) ((timeSteps - k) * (nonPeriodicEndAgent - nonPeriodicStartAgent + 1));
			avNegativeLocal /= (double) ((timeSteps - k) * (nonPeriodicEndAgent - nonPeriodicStartAgent + 1));
		}
		// Close off the individual observations
		try {
			if (computeMultiInfoCoherence) {
				miCalc.finaliseAddObservations();
			}
		} catch (Exception e) {
			// an exception would only be thrown if we changed the number of causal contributors here
			//  which simply will not happen. Just in case it does, we'll throw a runtime exception
			throw new RuntimeException("Number of causal contributors changed from intialisation to calculation!");
		}

		return localSep;
	}
	
	/**
	 * Computes local separable information for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method to be used for homogeneous agents only
	 *  since it assumes offsets are same for all destinations
	 *  and includes them all in the calculation
	 *  
	 * @param states multivariate time series 1st index is time, 2nd and 3rd index give the 2D agent number
	 * @param offsetOfDestFromSources offsets of the destination *from* causal information contributors.
	 * 		  (i.e. an offset of 1 means the destination is one index larger, or one to the right,
	 * 			than the source). 
	 *        sourcesOffsets is permitted to include 0, it will be ignored.
	 * @return multivariate time series of local separable information
	 *  values, indexed as per states.
	 */
	public double[][][] computeLocalFromPreviousObservations
		(int states[][][], int offsetOfDestFromSources[][]){
		
		return computeLocalFromPreviousObservations(states, offsetOfDestFromSources, false);
	}
	private double[][][] computeLocalFromPreviousObservations
		(int states[][][], int offsetOfDestFromSources[][], boolean cleanedOffsets){
		
		int[][] cleanedOffsetOfDestFromSources;
		if (cleanedOffsets) {
			cleanedOffsetOfDestFromSources = offsetOfDestFromSources;
		} else {
			cleanedOffsetOfDestFromSources = cleanOffsetOfDestFromSources(offsetOfDestFromSources);
			// This call made redundant by cleanOffsetSources:
			// confirmEnoughOffsetSources(othersOffsets);
		}

		int timeSteps = states.length;
		int numAgentRows = states[0].length;
		int numAgentColumns = states[0][0].length;

		int minRowOffset = MatrixUtils.min(cleanedOffsetOfDestFromSources, ROW_INDEX);
		int maxRowOffset = MatrixUtils.max(cleanedOffsetOfDestFromSources, ROW_INDEX);
		int nonPeriodicStartRow = Math.max(0, maxRowOffset);
		int nonPeriodicEndRow = numAgentRows - 1 + Math.min(0, minRowOffset);
		int minColumnOffset = MatrixUtils.min(cleanedOffsetOfDestFromSources, COLUMN_INDEX);
		int maxColumnOffset = MatrixUtils.max(cleanedOffsetOfDestFromSources, COLUMN_INDEX);
		int nonPeriodicStartColumn = Math.max(0, maxColumnOffset);
		int nonPeriodicEndColumn = numAgentColumns - 1 + Math.min(0, minColumnOffset);

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][][] localSep = new double[timeSteps][numAgentRows][numAgentColumns];
		average = 0.0;
		avPositiveLocal = 0.0;
		avNegativeLocal = 0.0;
		max = 0.0;
		min = 0.0;

		// Initialise and store the current previous value for each column
		int[][] pastVal = new int[numAgentRows][numAgentColumns];
		for (int r = 0; r < numAgentRows; r++) {
			for (int c = 0; c < numAgentColumns; c++) {
				pastVal[r][c] = 0;
				for (int p = 0; p < k; p++) {
					pastVal[r][c] *= base;
					pastVal[r][c] += states[p][r][c];
				}
			}
		}
		
		// Make a vector of the active and TE values for the coherence computation
		if (computeMultiInfoCoherence) {
			miCalc.startAddObservations();
		}
		double[] localActAndTes = new double[numSources + 1];
		
		int destVal, sourceVal;
		double logTerm;
		for (int t = k; t < timeSteps; t++) {
			for (int r = periodicBoundaryConditions ? 0 : nonPeriodicStartRow;
			 r < (periodicBoundaryConditions ? numAgentRows : nonPeriodicEndRow + 1);
			 r++) {
				for (int c = periodicBoundaryConditions ? 0 : nonPeriodicStartColumn;
				c < (periodicBoundaryConditions ? numAgentColumns : nonPeriodicEndColumn + 1);
				c++) {
					destVal = states[t][r][c];
					// Add in the local active information storage:
					logTerm = ( (double) nextPastCount[destVal][pastVal[r][c]] ) /
			  		  ( (double) nextCount[destVal] *
			  			(double) pastCount[pastVal[r][c]] );
					// Now account for the fact that we've just used counts rather than probabilities,
					//  and we've got two counts on the bottom but one count on the top:
					logTerm *= (double) observations;
					if (computeMultiInfoCoherence) {
						// Keep the local active info.
						// Adding natural logs, since we're going to normalize these anyway
						localActAndTes[0] = Math.log(logTerm);
					}
					for (int sIndex = 0; sIndex < numSources; sIndex++) {
						// Can safely do mod operations here - if not periodic boundary conditions
						//  and this source would have bounced over boundary, we would have skipped
						//  over this destination earlier
						sourceVal = states[t-1]
						        [(r-cleanedOffsetOfDestFromSources[sIndex][ROW_INDEX]+numAgentRows) % numAgentRows]
						        [(c-cleanedOffsetOfDestFromSources[sIndex][COLUMN_INDEX]+numAgentColumns) % numAgentColumns];
						localActAndTes[sIndex+1] = ((double) sourceNumValueNextPastCount[sIndex][sourceVal][destVal][pastVal[r][c]] / (double) sourceNumValuePastCount[sIndex][sourceVal][pastVal[r][c]]) /
				 			((double) nextPastCount[destVal][pastVal[r][c]] / (double) pastCount[pastVal[r][c]]);
						// Add in the local transfer entropy for each source:
						logTerm *= localActAndTes[sIndex+1];
						if (computeMultiInfoCoherence) {
							// Keep this local TE
							localActAndTes[sIndex+1] = Math.log(localActAndTes[sIndex+1]);
						}
					}
					// Add in these local active and TE values for the coherence calculation
					if (computeMultiInfoCoherence) {
						miCalc.addObservation(localActAndTes);
					}
					localSep[t][r][c] = Math.log(logTerm) / log_2;
					average += localSep[t][r][c];
					if (localSep[t][r][c] > 0.0) {
						avPositiveLocal += localSep[t][r][c];
					} else {
						avNegativeLocal += localSep[t][r][c];
					}
					if (localSep[t][r][c] > max) {
						max = localSep[t][r][c];
					} else if (localSep[t][r][c] < min) {
						min = localSep[t][r][c];
					}
					// Update the previous value:
					pastVal[r][c] -= maxShiftedValue[states[t-k][r][c]];
					pastVal[r][c] *= base;
					pastVal[r][c] += states[t][r][c];
				}
			}
		}		

		if (periodicBoundaryConditions) {
			average /= (double) ((timeSteps - k) * numAgentRows * numAgentColumns);
			avPositiveLocal /= (double) ((timeSteps - k) * numAgentRows * numAgentColumns);
			avNegativeLocal /= (double) ((timeSteps - k) * numAgentRows * numAgentColumns);
		} else {
			average /= (double) ((timeSteps - k) * (nonPeriodicEndRow - nonPeriodicStartRow + 1) *
					(nonPeriodicEndColumn - nonPeriodicStartColumn + 1));
			avPositiveLocal /= (double) ((timeSteps - k) * (nonPeriodicEndRow - nonPeriodicStartRow + 1) *
					(nonPeriodicEndColumn - nonPeriodicStartColumn + 1));
			avNegativeLocal /= (double) ((timeSteps - k) * (nonPeriodicEndRow - nonPeriodicStartRow + 1) *
					(nonPeriodicEndColumn - nonPeriodicStartColumn + 1));
		}
		// Close off the individual observations
		try {
			if (computeMultiInfoCoherence) {
				miCalc.finaliseAddObservations();
			}
		} catch (Exception e) {
			// an exception would only be thrown if we changed the number of causal contributors here
			//  which simply will not happen. Just in case it does, we'll throw a runtime exception
			throw new RuntimeException("Number of causal contributors changed from intialisation to calculation!");
		}

		return localSep;
	}

	/**
	 * Computes local separable information for the given
	 *  destination in the multivariate time series,
	 *  using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method is suitable for heterogeneous agents
	 *  
	 * @param states multivariate time series 1st index is time, 2nd index is agent number
	 * @param destCol index for the destination agent
	 * @param sourcesAbsolute indices for the source agents
	 * @return multivariate time series of local separable information
	 *  values for the given destination variable
	 */
	public double[] computeLocalFromPreviousObservations
		(int states[][], int destCol, int[] sourcesAbsolute){
		
		return computeLocalFromPreviousObservations(states, destCol, sourcesAbsolute, false);
	}
	private double[] computeLocalFromPreviousObservations
		(int states[][], int destCol, int[] sourcesAbsolute, boolean cleanedSources){

		int[] cleanedSourcesAbsolute;
		if (cleanedSources) {
			cleanedSourcesAbsolute = sourcesAbsolute;
		} else {
			cleanedSourcesAbsolute = cleanAbsoluteSources(sourcesAbsolute, destCol);
			// This call made redundant by cleanAbsoluteSources:
			// confirmEnoughAbsoluteSources(sourcesAbsolute, destCol);
		}
		
		int rows = states.length;
		// int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localSep = new double[rows];
		average = 0.0;
		avPositiveLocal = 0.0;
		avNegativeLocal = 0.0;
		max = 0.0;
		min = 0.0;

		// Initialise and store the current previous value for each column
		int pastVal = 0; 
		pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += states[p][destCol];
		}
		
		// Make a vector of the active and TE values for the coherence computation
		if (computeMultiInfoCoherence) {
			miCalc.startAddObservations();
		}
		double[] localActAndTes = new double[numSources + 1];

		int destVal, sourceVal;
		double logTerm;
		for (int r = k; r < rows; r++) {
			destVal = states[r][destCol];
			// Add in the active information storage
			logTerm = ( (double) nextPastCount[destVal][pastVal] ) /
	  		  ( (double) nextCount[destVal] *
	  			(double) pastCount[pastVal] );
			// Now account for the fact that we've just used counts rather than probabilities,
			//  and we've got two counts on the bottom but one count on the top:
			logTerm *= (double) observations;
			if (computeMultiInfoCoherence) {
				// Keep the local active info.
				// Adding natural logs, since we're going to normalize these anyway
				localActAndTes[0] = Math.log(logTerm);
			}
			for (int sIndex = 0; sIndex < numSources; sIndex++) {
				sourceVal = states[r-1][cleanedSourcesAbsolute[sIndex]];
				localActAndTes[sIndex+1] = ((double) sourceNumValueNextPastCount[sIndex][sourceVal][destVal][pastVal] / (double) sourceNumValuePastCount[sIndex][sourceVal][pastVal]) /
		 			((double) nextPastCount[destVal][pastVal] / (double) pastCount[pastVal]);
				// Add in the local transfer entropy for each source:
				logTerm *= localActAndTes[sIndex+1];
				if (computeMultiInfoCoherence) {
					// Keep this local TE
					localActAndTes[sIndex+1] = Math.log(localActAndTes[sIndex+1]);
				}
			}
			// Add in these local active and TE values for the coherence calculation
			if (computeMultiInfoCoherence) {
				miCalc.addObservation(localActAndTes);
			}
			localSep[r] = Math.log(logTerm) / log_2;
			average += localSep[r];
			if (localSep[r] > 0.0) {
				avPositiveLocal += localSep[r];
			} else {
				avNegativeLocal += localSep[r];
			}
			if (localSep[r] > max) {
				max = localSep[r];
			} else if (localSep[r] < min) {
				min = localSep[r];
			}
			// Update the previous value:
			pastVal -= maxShiftedValue[states[r-k][destCol]];
			pastVal *= base;
			pastVal += states[r][destCol];
		}

		average /= (double) (rows - k);
		avPositiveLocal /= (double) (rows - k);
		avNegativeLocal /= (double) (rows - k);
		
		// Close off the individual observations
		try {
			if (computeMultiInfoCoherence) {
				miCalc.finaliseAddObservations();
			}
		} catch (Exception e) {
			// an exception would only be thrown if we changed the number of causal contributors here
			//  which simply will not happen. Just in case it does, we'll throw a runtime exception
			throw new RuntimeException("Number of causal contributors changed from intialisation to calculation!");
		}

		return localSep;
	}

	/**
	 * Computes local separable information for the given
	 *  destination variable in a multivariate time series,
	 *  using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method is suitable for heterogeneous agents
	 *  
	 * @param states multivariate time series: 1st index is time,
	 * 	2nd and 3rd index give the 2D agent number
	 * @param destAgentRow the destination index
	 * @param destAgentColumn the destination index
	 * @param sourcesAbsolute array of the source indices
	 * @return multivariate time series of local separable information
	 *  values for the given destination
	 */
	public double[] computeLocalFromPreviousObservations
		(int states[][][], int destAgentRow, int destAgentColumn, int[][] sourcesAbsolute){
		
		return computeLocalFromPreviousObservations(states, destAgentRow, destAgentColumn, sourcesAbsolute, false);
	}
	private double[] computeLocalFromPreviousObservations
		(int states[][][], int destAgentRow, int destAgentColumn, int[][] sourcesAbsolute, boolean cleanedSources){

		int[][] cleanedSourcesAbsolute;
		if (cleanedSources) {
			cleanedSourcesAbsolute = sourcesAbsolute;
		} else {
			cleanedSourcesAbsolute = cleanAbsoluteSources(sourcesAbsolute, destAgentRow, destAgentColumn);
			// This call made redundant by cleanAbsoluteSources:
			// confirmEnoughAbsoluteSources(sourcesAbsolute, destCol);
		}
		
		int timeSteps = states.length;
		// int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localSep = new double[timeSteps];
		average = 0.0;
		avPositiveLocal = 0.0;
		avNegativeLocal = 0.0;
		max = 0.0;
		min = 0.0;

		// Initialise and store the current previous value for each column
		int pastVal = 0; 
		pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += states[p][destAgentRow][destAgentColumn];
		}
		
		// Make a vector of the active and TE values for the coherence computation
		if (computeMultiInfoCoherence) {
			miCalc.startAddObservations();
		}
		double[] localActAndTes = new double[numSources + 1];

		int destVal, sourceVal;
		double logTerm;
		for (int t = k; t < timeSteps; t++) {
			destVal = states[t][destAgentRow][destAgentColumn];
			// Add in the active information storage
			logTerm = ( (double) nextPastCount[destVal][pastVal] ) /
	  		  ( (double) nextCount[destVal] *
	  			(double) pastCount[pastVal] );
			// Now account for the fact that we've just used counts rather than probabilities,
			//  and we've got two counts on the bottom but one count on the top:
			logTerm *= (double) observations;
			if (computeMultiInfoCoherence) {
				// Keep the local active info.
				// Adding natural logs, since we're going to normalize these anyway
				localActAndTes[0] = Math.log(logTerm);
			}
			for (int sIndex = 0; sIndex < numSources; sIndex++) {
				sourceVal = states[t-1][cleanedSourcesAbsolute[sIndex][ROW_INDEX]]
				                       [cleanedSourcesAbsolute[sIndex][COLUMN_INDEX]];
				localActAndTes[sIndex+1] = ((double) sourceNumValueNextPastCount[sIndex][sourceVal][destVal][pastVal] / (double) sourceNumValuePastCount[sIndex][sourceVal][pastVal]) /
		 			((double) nextPastCount[destVal][pastVal] / (double) pastCount[pastVal]);
				// Add in the local transfer entropy for each source:
				logTerm *= localActAndTes[sIndex+1];
				if (computeMultiInfoCoherence) {
					// Keep this local TE
					localActAndTes[sIndex+1] = Math.log(localActAndTes[sIndex+1]);
				}
			}
			// Add in these local active and TE values for the coherence calculation
			if (computeMultiInfoCoherence) {
				miCalc.addObservation(localActAndTes);
			}
			localSep[t] = Math.log(logTerm) / log_2;
			average += localSep[t];
			if (localSep[t] > 0.0) {
				avPositiveLocal += localSep[t];
			} else {
				avNegativeLocal += localSep[t];
			}
			if (localSep[t] > max) {
				max = localSep[t];
			} else if (localSep[t] < min) {
				min = localSep[t];
			}
			// Update the previous value:
			pastVal -= maxShiftedValue[states[t-k][destAgentRow][destAgentColumn]];
			pastVal *= base;
			pastVal += states[t][destAgentRow][destAgentColumn];
		}

		average /= (double) (timeSteps - k);
		avPositiveLocal /= (double) (timeSteps - k);
		avNegativeLocal /= (double) (timeSteps - k);
		
		// Close off the individual observations
		try {
			if (computeMultiInfoCoherence) {
				miCalc.finaliseAddObservations();
			}
		} catch (Exception e) {
			// an exception would only be thrown if we changed the number of causal contributors here
			//  which simply will not happen. Just in case it does, we'll throw a runtime exception
			throw new RuntimeException("Number of causal contributors changed from intialisation to calculation!");
		}

		return localSep;
	}

	/**
	 * Standalone routine to 
	 * compute local separable info across a 2D multivariate time series
	 *  of the states of homogeneous agents
	 * Return a 2D spatiotemporal multivariate time series of local values.
	 * First history rows are zeros.
	 * This method to be called for homogeneous agents only
	 * since it assumes all destinations have the same source offsets
	 * and includes them all in the calculation.
	 * 
	 * @param states multivariate time series: first index is time,
	 *  second is variable index.
	 * @param offsetOfDestFromSources offsets of the destination *from* causal information contributors.
	 * 		  (i.e. an offset of 1 means the destination is one index larger, or one to the right,
	 * 			than the source). 
	 *        sourcesOffsets is permitted to include 0, it will be ignored.
	 * @return multivariate time series of local separable information,
	 *  indexed as per states.
	 */
	public double[][] computeLocal(int states[][], int[] offsetOfDestFromSources) {
		
		initialise();
		int[] cleanedSourcesOffsets = cleanOffsetOfDestFromSources(offsetOfDestFromSources);
		addObservations(states, cleanedSourcesOffsets, true);
		return computeLocalFromPreviousObservations(states, cleanedSourcesOffsets, true);
	}

	/**
	 * Standalone routine to 
	 * compute local separable info across a 3D multivariate time series
	 *  of the states of homogeneous agents.
	 * Return a 3D spatiotemporal array of local values.
	 * First history rows are zeros.
	 * This method to be called for homogeneous agents only
	 * since it assumes all destinations have the same source offsets
	 * and includes them all in the calculation.
	 * 
	 * @param states multivariate time series: 1st index is time, 2nd and 3rd are agent indices
	 * @param offsetOfDestFromSources offsets of the destination *from* causal information contributors.
	 * 		  (i.e. an offset of 1 means the destination is one index larger, or one to the right,
	 * 			than the source). 
	 *        sourcesOffsets is permitted to include 0, it will be ignored.
	 * @return multivariate time series of local separable information,
	 *  indexed as per states.
	 */
	public double[][][] computeLocal(int states[][][], int[][] offsetOfDestFromSources) {
		
		initialise();
		int[][] cleanedSourcesOffsets = cleanOffsetOfDestFromSources(offsetOfDestFromSources);
		addObservations(states, cleanedSourcesOffsets, true);
		return computeLocalFromPreviousObservations(states, cleanedSourcesOffsets, true);
	}

	/**
	 * Standalone routine to 
	 * compute average separable information across a 2D multivariate time series
	 *  of the states of homogeneous agents.
	 * Return the average.
	 * This method to be called for homogeneous agents only
	 * since it assumes all destinations have the same source offsets
	 * and includes them all in the calculation.
	 * 
	 * @param states multivariate time series: first index is time,
	 *  second is variable number
	 * @param sourceOffsets column offsets for causal info contributors
	 *  from each destination
	 * @return average separable information for these
	 *  observations.
	 */
	public double computeAverageLocal(int states[][], int[] sourceOffsets) {
		
		initialise();
		addObservations(states, sourceOffsets);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute average separable information across a 3D multivariate time series
	 *  of the states of homogeneous agents.
	 * Return the average separable info.
	 * This method to be called for homogeneous agents only
	 * since it assumes all destinations have the same source offsets
	 * and includes them all in the calculation.
	 * 
	 * @param states multivariate time series: 1st index is time, 2nd and 3rd are agent indices
	 * @param sourceOffsets agent offsets for causal info contributors. 1st index points to 
	 *  an array of two elements for the row and column offsets.
	 * @return average separable information for these
	 *  observations.
	 */
	public double computeAverageLocal(int states[][][], int[][] sourceOffsets) {
		
		initialise();
		addObservations(states, sourceOffsets);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute local separable information for a specific destination
	 * in a 2D multivariate time series
	 *  of the states.
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros.
	 * This method suitable for heterogeneous agents
	 * 
	 * @param states multivariate time series: first index is time,
	 *  second is variable number
	 * @param destCol column index for the destination agent to consider
	 * @param sourcesAbsolute column indices for causal info contributors
	 *  to this destination
	 * @return time series of local separable information
	 */
	public double[] computeLocal(int states[][], int destCol, int[] sourcesAbsolute) {
		
		initialise();
		int[] cleanedSources = cleanAbsoluteSources(sourcesAbsolute, destCol);
		addObservations(states, destCol, cleanedSources, true);
		return computeLocalFromPreviousObservations(states, destCol, cleanedSources, true);
	}

	/**
	 * Standalone routine to 
	 * compute local separable information for a specific destination
	 * in a 3D multivariate time series
	 *  of the states.
	 * First history rows are zeros.
	 * This method suitable for heterogeneous agents.
	 * 
	 * @param states multivariate time series: 1st index is time, 2nd and 3rd are agent indices
	 * @param destAgentRow row index for the destination agent
	 * @param destAgnetColumn column index for the destination agent
	 * @param sourcesAbsolute absolute indices for causal info contributors to this destination
	 * @return time series of local separable information
	 */
	public double[] computeLocal(int states[][][], int destAgentRow,
			int destAgentColumn, int[][] sourcesAbsolute) {
		
		initialise();
		int[][] cleanedSources = cleanAbsoluteSources(sourcesAbsolute, destAgentRow, destAgentColumn);
		addObservations(states, destAgentRow, destAgentColumn, cleanedSources, true);
		return computeLocalFromPreviousObservations(states, destAgentRow, destAgentColumn, cleanedSources, true);
	}

	/**
	 * Standalone routine to 
	 * compute average separable information for a specific destination
	 * in a 2D multivariate time series
	 *  of the states.
	 * Returns the average.
	 * This method suitable for heterogeneous agents
	 * 
	 * @param states multivariate time series: first index is time,
	 *  second is variable number
	 * @param destCol column index for the destination agent
	 * @param sourcesAbsolute column indices for causal info contributors
	 * @return average separable information for this destination
	 */
	public double computeAverageLocal(int states[][], int destCol, int[] sourcesAbsolute) {
		
		initialise();
		addObservations(states, destCol, sourcesAbsolute);
		return computeAverageLocalOfObservations();
	}
	
	/**
	 * Standalone routine to 
	 * compute average separable information for a specific destination
	 * in a 3D multivariate time series
	 *  of the states.
	 * Returns the average.
	 * This method suitable for heterogeneous agents
	 * 
	 * @param states multivariate time series: 1st index is time, 2nd and 3rd are agent indices
	 * @param destAgentRow row index for the destination agent
	 * @param destAgnetColumn column index for the destination agent
	 * @param sourcesAbsolute absolute indices for causal info contributors to this destination
	 * @return average separable information for this destination
	 */
	public double computeAverageLocal(int states[][][], int destAgentRow,
			int destAgentColumn, int[][] sourcesAbsolute) {
		
		initialise();
		addObservations(states, destAgentRow, destAgentColumn, sourcesAbsolute);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Returns the average of the positive components from the
	 *  last computation.
	 * The average is taken over all components, whether positive
	 *  or negative, such that:
	 *  getLastAveragePositive() + getLastAverageNegative()
	 *   == getLastAverage()
	 * 
	 * @return
	 */
	public double getLastAveragePositive() {
		return avPositiveLocal;
	}
	
	/**
	 * Returns the average of the negative components from the
	 *  last computation.
	 * The average is taken over all components, whether positive
	 *  or negative, such that:
	 *  getLastAveragePositive() + getLastAverageNegative()
	 *   == getLastAverage()
	 * 
	 * @return
	 */
	public double getLastAverageNegative() {
		return avNegativeLocal;
	}

	/**
	 * Counts the information contributors to this node which
	 * are not equal to the node itself (offset 0)
	 * 
	 * @param sourcesOffsets information contributors
	 * @return count of sourcesOffsets with any 0 entries removed from the array
	 */
	public static int countOfOffsetSources(int[] sourcesOffsets) {
		int countOfSources = 0;
		for (int index = 0; index < sourcesOffsets.length; index++) {
			if (sourcesOffsets[index] != 0) {
				countOfSources++;
			}
		}
		return countOfSources;
	}

	/**
	 * Counts the information contributors to this node which
	 * are not equal to the node itself (offset (0,0))
	 * 
	 * @param sourcesOffsets information contributors
	 * @return count of sourcesOffsets with any (0,0) entries removed from the 2D array
	 */
	public static int countOfOffsetSources(int[][] sourcesOffsets) {
		int countOfSources = 0;
		for (int index = 0; index < sourcesOffsets.length; index++) {
			if ((sourcesOffsets[index][ROW_INDEX] != 0) && (sourcesOffsets[index][COLUMN_INDEX] != 0)) {
				countOfSources++;
			}
		}
		return countOfSources;
	}

	/**
	 * Counts the information contributors to the dest which
	 * are not equal to the node itself 
	 * 
	 * @param sources array of source indices
	 * @param dest index of dest
	 * @return count of sources with any entries equal to dest removed.
	 */
	public static int countOfAbsoluteSources(int[] sources, int dest) {
		int countOfSources = 0;
		for (int index = 0; index < sources.length; index++) {
			if (sources[index] != dest) {
				countOfSources++;
			}
		}
		return countOfSources;
	}

	/**
	 * Counts the information contributors to the dest which
	 * are not equal to the node itself 
	 * 
	 * @param sources array of arrays of row and column indices
	 * @param destAgentRow row of dest variable
	 * @param destAgentColumn column of dest variable
	 * @return count of entries in sources not equal to the dest indices
	 */
	public static int countOfAbsoluteSources(int[][] sources, int destAgentRow, int destAgentColumn) {
		int countOfSources = 0;
		for (int index = 0; index < sources.length; index++) {
			if ((sources[index][ROW_INDEX] != destAgentRow) &&
					(sources[index][COLUMN_INDEX] != destAgentColumn)) {
				countOfSources++;
			}
		}
		return countOfSources;
	}

	/**
	 * Check that the supplied array of offsets as sources
	 *  is long enough compared to our expectation
	 * 
	 * @param sourcesOffsets array of source offsets
	 * @return whether it is as long as expected
	 * @throws RuntimeException if it is not
	 */
	public boolean confirmEnoughOffsetSources(int[] sourcesOffsets) {
		if (countOfOffsetSources(sourcesOffsets) != numSources) {
			throw new RuntimeException("Incorrect number of sources in offsets");
		}
		return true;
	}

	/**
	 * Check that the supplied array of offsets as sources
	 *  is long enough compared to our expectation
	 * 
	 * @param sourcesOffsets array of 2D source offsets
	 * @return whether it is as long as expected
	 * @throws RuntimeException if it is not
	 */
	public boolean confirmEnoughOffsetSources(int[][] sourcesOffsets) {
		if (countOfOffsetSources(sourcesOffsets) != numSources) {
			throw new RuntimeException("Incorrect number of sources in offsets");
		}
		return true;
	}

	/**
	 * Check that the supplied array of absolutes as sources
	 *  is long enough compared to our expectation
	 * 
	 * @param sourcesAbsolute array of source indices
	 * @param dest dest index
	 * @return whether it is as long as expected
	 * @throws RuntimeException if it is not
	 */
	public boolean confirmEnoughAbsoluteSources(int[] sourcesAbsolute, int dest) {
		if (countOfAbsoluteSources(sourcesAbsolute, dest) != numSources) {
			throw new RuntimeException("Incorrect number of sources in absolutes");
		}
		return true;
	}
	
	/**
	 * Check that the supplied array of absolutes as sources
	 *  is long enough compared to our expectation
	 * 
	 * @param sourcesAbsolute array of 2D source indices
	 * @param destAgentRow row of dest variable
	 * @param destAgentColumn column of dest variable
	 * @return whether it is as long as expected
	 * @throws RuntimeException if it is not
	 */
	public boolean confirmEnoughAbsoluteSources(int[][] sourcesAbsolute, int destAgentRow, int destAgentColumn) {
		if (countOfAbsoluteSources(sourcesAbsolute, destAgentRow, destAgentColumn) != numSources) {
			throw new RuntimeException("Incorrect number of sources in absolutes");
		}
		return true;
	}

	/**
	 * Returns the information contributors to this node which
	 * are not equal to the node itself (offset 0).
	 * Checks that there are enough sources.
	 * 
	 * @param sourcesOffsets array of source offsets
	 * @return sourcesOffsets with these entries removed
	 *  from the array
	 */
	public int[] cleanOffsetOfDestFromSources(int[] sourcesOffsets) {
		int[] cleaned = new int[numSources];
		int countOfSources = 0;
		for (int index = 0; index < sourcesOffsets.length; index++) {
			if (sourcesOffsets[index] != 0) {
				if (countOfSources == numSources) {
					// We've already taken all the sources we expected
					countOfSources++;
					break;
				}
				cleaned[countOfSources] = sourcesOffsets[index];
				countOfSources++;
			}
		}
		if (countOfSources < numSources) {
			throw new RuntimeException("Too few sources in offsets");
		} else if (countOfSources > numSources) {
			throw new RuntimeException("Too many sources in offsets");
		}
		return cleaned;
	}
	
	/**
	 * Returns the information contributors to this node which
	 * are not equal to the node itself (offset (0,0)).
	 * Checks that there are enough sources.
	 * 
	 * @param sourcesOffsets 2D source offsets. 1st dimension is source index, 2nd index is for
	 *   1st or 2nd index of source indice pair
	 * @return sourcesOffsets with these entries removed
	 *  from the array
	 */
	public int[][] cleanOffsetOfDestFromSources(int[][] sourcesOffsets) {
		int[][] cleaned = new int[numSources][2];
		int countOfSources = 0;
		for (int index = 0; index < sourcesOffsets.length; index++) {
			if ((sourcesOffsets[index][0] != 0) || (sourcesOffsets[index][1] != 0)){
				if (countOfSources == numSources) {
					// We've already taken all the sources we expected
					countOfSources++;
					break;
				}
				// copy both indices of the source
				cleaned[countOfSources][ROW_INDEX] = sourcesOffsets[index][ROW_INDEX];
				cleaned[countOfSources][COLUMN_INDEX] = sourcesOffsets[index][COLUMN_INDEX];
				countOfSources++;
			}
		}
		if (countOfSources < numSources) {
			throw new RuntimeException("Too few sources in offsets");
		} else if (countOfSources > numSources) {
			throw new RuntimeException("Too many sources in offsets");
		}
		return cleaned;
	}

	/**
	 * Returns the information contributors to the dest which
	 * are not equal to the node itself (offset 0).
	 * Checks that there are enough other information contributors.
	 * 
	 * @param sources array of source indices
	 * @param dest dest index
	 * @return sources with these entries removed
	 *  from the array
	 */
	public int[] cleanAbsoluteSources(int[] sources, int dest) {
		int[] cleaned = new int[numSources];
		int countOfSources = 0;
		for (int index = 0; index < sources.length; index++) {
			if (sources[index] != dest) {
				if (countOfSources == numSources) {
					// We've already taken all the other info 
					//  contributors we expected
					countOfSources++;
					break;
				}
				cleaned[countOfSources] = sources[index];
				countOfSources++;
			}
		}
		if (countOfSources < numSources) {
			throw new RuntimeException("Too few sources in absolutes");
		} else if (countOfSources > numSources) {
			throw new RuntimeException("Too many sources in absolutes");
		}
		return cleaned;
	}
	
	/**
	 * Returns the information contributors to the dest which
	 * are not equal to the node itself (offset 0).
	 * Checks that there are enough other information contributors.
	 * 
	 * @param sources 2D array of source indices
	 * @param dest dest index
	 * @return sources with these entries removed
	 *  from the array
	 */
	public int[][] cleanAbsoluteSources(int[][] sources, int destAgentRow, int destAgentColumn) {
		int[][] cleaned = new int[numSources][2];
		int countOfSources = 0;
		for (int index = 0; index < sources.length; index++) {
			if ((sources[index][ROW_INDEX] != destAgentRow) &&
					(sources[index][COLUMN_INDEX] != destAgentColumn)){
				if (countOfSources == numSources) {
					// We've already taken all the other info 
					//  contributors we expected
					countOfSources++;
					break;
				}
				cleaned[countOfSources][ROW_INDEX] = sources[index][ROW_INDEX];
				cleaned[countOfSources][COLUMN_INDEX] = sources[index][COLUMN_INDEX];
				countOfSources++;
			}
		}
		if (countOfSources < numSources) {
			throw new RuntimeException("Too few sources in absolutes");
		} else if (countOfSources > numSources) {
			throw new RuntimeException("Too many sources in absolutes");
		}
		return cleaned;
	}
	
	/**
	 * Whether we use periodic boundary conditions when considering
	 * spatiotemporal data for homogeneous variables
	 * 
	 * @return if this is the case
	 */
	public boolean isPeriodicBoundaryConditions() {
		return periodicBoundaryConditions;
	}
	/**
	 * Set whether to use periodic boundary conditions
	 * 
	 * @param periodicBoundaryConditions as above
	 */
	public void setPeriodicBoundaryConditions(boolean periodicBoundaryConditions) {
		this.periodicBoundaryConditions = periodicBoundaryConditions;
	}
	
	/**
	 * Whether the calculator should gather observations to
	 *  compute a multi-info between the local TEs and
	 *  AIS values.
	 * See Lizier et al, 2012 "Coherent information structure in complex computation"
	 *   
	 */
	public boolean isComputeMultiInfoCoherence() {
		return computeMultiInfoCoherence;
	}
	
	public void setComputeMultiInfoCoherence(String miCoherenceCalculatorClass,
			Properties props) throws Exception {
		computeMultiInfoCoherence = true;
		try {
			miCalc = (MultiInfoCalculator) Class.forName(miCoherenceCalculatorClass).newInstance();
		} catch (Exception e) {
			throw new RuntimeException("Cannot initiate class " + miCoherenceCalculatorClass +
					" as the MultiInfoCalculator class inside SeparableInfoCalculator");
		}
		miCalc.setDebug(debug);
		// Pass all the properties onto the calculator, and let it
		//  determine which ones are relevant
		for (Object propertyNameObj : props.keySet()) {
			String propertyName = (String) propertyNameObj;
			String propertyValue = props.getProperty(propertyName);
			miCalc.setProperty(propertyName, propertyValue);
		}
		if (miCalc == null) {
			throw new Exception("Calculator was not initialised with a multi-info calculator to compute the coherence");
		}
	}
	
	public void clearComputeMultiInfoCoherence() {
		computeMultiInfoCoherence = false;
	}
	
	public double computeMultiInfoCoherence() throws Exception {
		if (!computeMultiInfoCoherence) {
			throw new Exception("Calculator was not set to track coherence before separable info was calculated");
		}
		return miCalc.computeAverageLocalOfObservations();
	}
	
	/**
	 * 
	 * @return whether the calculator can compute the multi info
	 *  coherence of computation from the averageLocalOfObservations
	 *  method.
	 */
	public boolean canComputeMultiInfoCoherenceFromAverageOfObservations() {
		return true;
	}
	
	/**
	 * Allows reclaiming of some vital memory
	 */
	public void resetMultiInfoCoherenceCalculator() {
		miCalc.initialise(numSources + 1);
	}
	
	@Override
	public void setDebug(boolean debug) {
		if (computeMultiInfoCoherence) {
			miCalc.setDebug(debug);
		}
		super.setDebug(debug);
	}

}
