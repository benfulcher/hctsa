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

import infodynamics.utils.AnalyticMeasurementDistribution;
import infodynamics.utils.AnalyticNullDistributionComputer;
import infodynamics.utils.ChiSquareMeasurementDistribution;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

/**
 * <p>Implements <b>transfer entropy</b>
 * for univariate discrete time-series data.
 * That is, it is applied to <code>int[]</code> data, indexed
 * by time.
 * See Schreiber below for the definition of transfer entropy,
 * and Lizier et al. for the definition of local transfer entropy.
 * Specifically, this class implements the pairwise or <i>apparent</i>
 * transfer entropy; i.e. we compute the transfer that appears to
 * come from a single source variable, without examining any other
 * potential sources
 * (see Lizier et al, PRE, 2008).</p>
 *  
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * <ol>
 * 		<li>Construct the calculator via {@link #TransferEntropyCalculatorDiscrete(int, int)}
 * 			or {@link #TransferEntropyCalculatorDiscrete(int, int, int)};</li>
 *		<li>Initialise the calculator using
 *			{@link #initialise()};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using one or more calls to
 * 			the set of {@link #addObservations(int[], int[])} methods, then</li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average TE: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>the local TE values for these samples: {@link #computeLocalOfPreviousObservations()}</li>
 * 				<li>local TE values for a specific set of samples: e.g.
 * 				{@link #computeLocalFromPreviousObservations(int[], int[])} etc.</li>
 * 				<li>the distribution of TE values under the null hypothesis
 * 					of no relationship between source and
 * 					destination values: {@link #computeSignificance(int)} or
 * 					{@link #computeSignificance(int[][])}.</li>
 * 			</ul>
 * 		</li>
 * 		<li>As an alternative to steps 3 and 4, the user may undertake
 * 			standalone computation from a single set of observations, via
 *  		e.g.: {@link #computeLocal(int[], int[])} or
 *  		{@link #computeAverageLocal(int[][], int)}.</li>
 * 		<li>
 * 		Return to step 2 to re-use the calculator on a new data set.
 * 		</li>
 * 	</ol>
 * 
 * TODO Add arbitrary source-dest delay
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. Schreiber, <a href="http://dx.doi.org/10.1103/PhysRevLett.85.461">
 * "Measuring information transfer"</a>,
 *  Physical Review Letters 85 (2) pp.461-464, 2000.</li>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href="http://dx.doi.org/10.1103/PhysRevE.77.026110">
 *  "Local information transfer as a spatiotemporal filter for complex systems"</a>
 *  Physical Review E 77, 026110, 2008.</li>
 * </ul>
 * 
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 */
public class TransferEntropyCalculatorDiscrete extends ContextOfPastMeasureCalculatorDiscrete 
	implements ChannelCalculatorDiscrete, AnalyticNullDistributionComputer {

	/**
	 * Counts of (source,dest_next,dest_embedded_past) tuples
	 */
	protected int[][][] sourceNextPastCount = null;	// count for (source[n],dest[n+1],dest[n]^k) tuples
	/**
	 * Counts of (source,dest_embedded_past) tuples
	 */
	protected int[][] sourcePastCount = null;			// count for (source[n],dest[n]^k) tuples
	/**
	 * Whether to assume periodic boundary conditions for channels across
	 *  the boundary of the multidimensional
	 *  calls supplying observations, e.g.
	 *  {@link #addObservations(int[][], int)} calls
	 */
	protected boolean periodicBoundaryConditions = true;

	/**
	 * Embedding length of the source variable.
	 * This is "l" in Schreiber's notation.
	 */
	protected int sourceHistoryEmbedLength = 1;
	
	/**
	 * A cached value of base^sourceHistoryEmbedLength
	 */
	protected int base_power_l = 1;
	
	/**
	 * A cached value of each discrete value left shifted (in "base" counting) by (sourceHistoryEmbedLength-1).
	 */
	protected int[] maxShiftedSourceValue = null; // states * (base^(sourceHistoryEmbedLength-1))

	/**
	 * First time step at which we can take an observation
	 *  (needs to account for an embedding in the previous steps)
	 */
	protected int startObservationTime = 1;
	
	/**
	 * Tracks whether the measure has been computed since the last initialisation
	 */
	protected boolean estimateComputed = false;
	
	/**
	 * User was formerly forced to create new instances through this factory method.
	 * Retained for backwards compatibility.
	 * 
	 * @param base
	 * @param destHistoryEmbedLength
	 * 
	 * @return a new TransferEntropyCalculator object
	 * @deprecated
	 */
	public static TransferEntropyCalculatorDiscrete newInstance(int base, int destHistoryEmbedLength) {
		
		return new TransferEntropyCalculatorDiscrete(base, destHistoryEmbedLength);

		// Old code for an attempted optimisation:
		/*
		if (isPowerOf2(base)) {
			return new ApparentTransferEntropyCalculatorBase2(base, history);
		} else {
			return new ApparentTransferEntropyCalculator(base, history);
		}
		*/
	}
	
	/**
	 * Create a new TE calculator for the given base and destination history embedding length.
	 * 
	 * @param base number of symbols for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param destHistoryEmbedLength embedded history length of the destination to condition on -
	 *        this is k in Schreiber's notation.
	 */
	public TransferEntropyCalculatorDiscrete(int base, int destHistoryEmbedLength) {

		super(base, destHistoryEmbedLength);
		base_power_l = MathsUtils.power(base, sourceHistoryEmbedLength);
		
		// Create constants for tracking sourceValues
		maxShiftedSourceValue = new int[base];
		for (int v = 0; v < base; v++) {
			maxShiftedSourceValue[v] = v;
		}

		// Create storage for extra counts of observations
		sourceNextPastCount = new int[base_power_l][base][base_power_k];
		sourcePastCount = new int[base_power_l][base_power_k];
		
		// Which time step do we start taking observations from?
		// Normally this is k (to allow k previous time steps)
		//  but if k==0 (becoming a lagged MI), it's 1.
		startObservationTime = Math.max(k, 1);
	}

	/**
	 * Create a new TE calculator for the given base, destination and source history embedding lengths.
	 * 
	 * @param base number of quantisation levels for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param destHistoryEmbedLength embedded history length of the destination to condition on -
	 *        this is k in Schreiber's notation.
	 * @param sourceHistoryEmbeddingLength embedded history length of the source to include -
	 *        this is l in Schreiber's notation.
	 */
	public TransferEntropyCalculatorDiscrete(int base, int destHistoryEmbedLength, int sourceHistoryEmbeddingLength) {

		super(base, destHistoryEmbedLength);
		this.sourceHistoryEmbedLength = sourceHistoryEmbeddingLength;
		base_power_l = MathsUtils.power(base, sourceHistoryEmbedLength);
		
		// Check that we can convert the history value into an integer ok: 
		if (sourceHistoryEmbedLength > Math.log(Integer.MAX_VALUE) / log_base) {
			throw new RuntimeException("Base and source history combination too large");
		}

		// Create constants for tracking sourceValues
		maxShiftedSourceValue = new int[base];
		for (int v = 0; v < base; v++) {
			maxShiftedSourceValue[v] = v * MathsUtils.power(base, sourceHistoryEmbedLength-1);
		}

		// Create storage for extra counts of observations
		sourceNextPastCount = new int[base_power_l][base][base_power_k];
		sourcePastCount = new int[base_power_l][base_power_k];
		
		// Which time step do we start taking observations from?
		// Normally this is k (to allow k previous time steps)
		//  but if k==0 (becoming a lagged MI), it's 1.
		// We also allow for source embeddings here too
		startObservationTime = Math.max(Math.max(k, sourceHistoryEmbedLength), 1);
	}

	@Override
	public void initialise(){
		super.initialise();
		estimateComputed = false;
		
		MatrixUtils.fill(sourceNextPastCount, 0);
		MatrixUtils.fill(sourcePastCount, 0);
	}
	
	@Override
	public void addObservations(int[] source, int[] dest) {
		int rows = dest.length;
		// increment the count of observations:
		observations += (rows - startObservationTime); 
		
		// Initialise and store the current previous value
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += dest[startObservationTime - k + p];
		}
		int sourceVal = 0;
		for (int p = 0; p < sourceHistoryEmbedLength; p++) {
			sourceVal *= base;
			sourceVal += source[startObservationTime - sourceHistoryEmbedLength + p];
		}
		
		// 1. Count the tuples observed
		int destVal;
		for (int r = startObservationTime; r < rows; r++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			destVal = dest[r];
			sourceNextPastCount[sourceVal][destVal][pastVal]++;
			sourcePastCount[sourceVal][pastVal]++;
			nextPastCount[destVal][pastVal]++;
			pastCount[pastVal]++;
			nextCount[destVal]++;
			// Update the previous values:
			if (k > 0) {
				pastVal -= maxShiftedValue[dest[r-k]];
				pastVal *= base;
				pastVal += dest[r];
			}
			sourceVal -= maxShiftedSourceValue[source[r-sourceHistoryEmbedLength]];
			sourceVal *= base;
			sourceVal += source[r];
		}
	}

	/**
 	 * Add observations for a single source-destination pair 
 	 *  to our estimates of the pdfs.
 	 *  
	 * @param source source time-series
	 * @param dest destination time-series. 
	 *  Must be same length as source
	 * @param valid time-series of whether the signals
	 *  at the given time should be considered valid 
	 *  and added to our PDFs
	 */
	public void addObservations(int[] source, int[] dest, boolean[] valid) {
		int rows = dest.length;
		
		// Initialise and store the current previous value
		int pastVal = 0;
		// We can take an observation if timeSinceLastDestInvalid > k
		//  and timeSinceLastSourceInvalid > sourceHistoryEmbedLength.
		// In preparation for introducing a source-dest delay,
		//  we're using two different time variables here.
		int timeSinceLastDestInvalid = k;
		int timeSinceLastSourceInvalid = sourceHistoryEmbedLength;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += dest[startObservationTime - k + p];
			if (!valid[startObservationTime - k + p]) {
				// Time from startObservationTime backwards to this step
				timeSinceLastDestInvalid =  k - p - 1;
			}
		}
		int sourceVal = 0;
		for (int p = 0; p < sourceHistoryEmbedLength; p++) {
			sourceVal *= base;
			sourceVal += source[startObservationTime - sourceHistoryEmbedLength + p];
			if (!valid[startObservationTime - sourceHistoryEmbedLength + p]) {
				// Time from startObservationTime backwards to this step
				timeSinceLastSourceInvalid =  sourceHistoryEmbedLength - p - 1;
			}
		}
		
		// 1. Count the tuples observed
		int destVal;
		for (int r = startObservationTime; r < rows; r++) {
			timeSinceLastDestInvalid++;
			timeSinceLastSourceInvalid++;
			// Pre-condition now:
			//  timeSinceLastDestInvalid holds the time from r back to
			//   the last valid destination point (not including current destination)
			//  timeSinceLastSourceInvalid holds the time from r back to
			//   the last valid source point (not including new source value
			
			if (!valid[r]) {
				timeSinceLastDestInvalid = 0;
			} else if ((timeSinceLastDestInvalid > k) &&
					    (timeSinceLastSourceInvalid > sourceHistoryEmbedLength)) {
				// Add to the count for this particular transition:
				// (cell's assigned as above)
				destVal = dest[r];
				sourceNextPastCount[sourceVal][destVal][pastVal]++;
				sourcePastCount[sourceVal][pastVal]++;
				nextPastCount[destVal][pastVal]++;
				pastCount[pastVal]++;
				nextCount[destVal]++;
				observations++;
			}
			// Update the previous values:
			if (k > 0) {
				pastVal -= maxShiftedValue[dest[r-k]];
				pastVal *= base;
				pastVal += dest[r];
			}
			sourceVal -= maxShiftedSourceValue[source[r-sourceHistoryEmbedLength]];
			sourceVal *= base;
			sourceVal += source[r];
			if (!valid[r]) {
				timeSinceLastSourceInvalid = 0;
			}
		}
	}

	/**
 	 * Add observations for a single source-destination pair
 	 *  to our estimates of the pdfs.
	 * Start and end time are the (inclusive) indices within which to add the observations.
	 * The start time is from the earliest of the k historical values of the destination (inclusive),
	 *  the end time is the last destination time point to add in.
	 *  
	 * @param source source time-series
	 * @param dest destination time-series. 
	 *  Must be same length as source
	 * @param startTime earliest time that we may extract embedded history from
	 * @param endTime last destination (next) time point to add in
	 * 
	 */
	public void addObservations(int[] source, int[] dest, int startTime, int endTime) {
		// increment the count of observations:
		observations += (endTime - startTime) - startObservationTime + 1; 
		
		// Initialise and store the current previous value
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += dest[startTime + startObservationTime - k + p];
		}
		int sourceVal = 0;
		for (int p = 0; p < sourceHistoryEmbedLength; p++) {
			sourceVal *= base;
			sourceVal += source[startTime + startObservationTime - sourceHistoryEmbedLength + p];
		}
		
		// 1. Count the tuples observed
		int destVal;
		for (int r = startTime + startObservationTime; r <= endTime; r++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			destVal = dest[r];
			sourceNextPastCount[sourceVal][destVal][pastVal]++;
			sourcePastCount[sourceVal][pastVal]++;
			nextPastCount[destVal][pastVal]++;
			pastCount[pastVal]++;
			nextCount[destVal]++;
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[dest[r-k]];
				pastVal *= base;
				pastVal += dest[r];
			}
			sourceVal -= maxShiftedSourceValue[source[r-sourceHistoryEmbedLength]];
			sourceVal *= base;
			sourceVal += source[r];
		}
	}

	/**
 	 * Add observations in to our estimates of the PDFs,
 	 * from a multivariate time-series.
 	 * This call suitable only for homogeneous agents, as all
 	 *  variable pairs separated by j column will contribute to the PDFs.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. source is column i-j, dest is column i: we
	 *  compute transfer is j cells to the right, using observations
	 *  across all column pairs separated by j) 
	 */
	public void addObservations(int states[][], int j) {
		int timeSteps = states.length;
		int agents = states[0].length;
		// increment the count of observations:
		if (periodicBoundaryConditions) {
			observations += (timeSteps - startObservationTime)*agents; 			
		} else {
			observations += (timeSteps - startObservationTime)*(agents - Math.abs(j));
		}
		
		// Initialise and store the current previous and source value for each column
		int[] pastVal = new int[agents];
		for (int c = 0; c < agents; c++) {
			pastVal[c] = 0;
			for (int p = 0; p < k; p++) {
				pastVal[c] *= base;
				pastVal[c] += states[startObservationTime - k + p][c];
			}
		}
		int[] sourceVal = new int[agents];
		for (int c = 0; c < agents; c++) {
			sourceVal[c] = 0;
			int sourceAgent = c-j;
			if ((sourceAgent < 0) || (sourceAgent >= agents)) {
				// Source agent is out of bounds unless we are using periodic boundary conditions
				if (periodicBoundaryConditions) {
					sourceAgent = (sourceAgent+agents) % agents;
				} else {
					// Don't add this to our observations
					continue;
				}
			}
			for (int p = 0; p < sourceHistoryEmbedLength; p++) {
				sourceVal[c] *= base;
				sourceVal[c] += states[startObservationTime - sourceHistoryEmbedLength + p][sourceAgent];
			}
		}
		
		// 1. Count the tuples observed
		int destVal;
		for (int r = startObservationTime; r < timeSteps; r++) {
			for (int c = 0; c < agents; c++) {
				// Add to the count for this particular transition:
				// (cell's assigned as above)
				int sourceAgent = c-j;
				if ((sourceAgent < 0) || (sourceAgent >= agents)) {
					// Source agent is out of bounds unless we are using periodic boundary conditions
					if (periodicBoundaryConditions) {
						sourceAgent = (sourceAgent+agents) % agents;
					} else {
						// Don't add this to our observations
						continue;
					}
				}
				destVal = states[r][c];
				sourceNextPastCount[sourceVal[c]][destVal][pastVal[c]]++;
				sourcePastCount[sourceVal[c]][pastVal[c]]++;
				nextPastCount[destVal][pastVal[c]]++;
				pastCount[pastVal[c]]++;
				nextCount[destVal]++;
				// Update the previous value:
				if (k > 0) {
					pastVal[c] -= maxShiftedValue[states[r-k][c]];
					pastVal[c] *= base;
					pastVal[c] += states[r][c];
				}
				sourceVal[c] -= maxShiftedSourceValue[states[r-sourceHistoryEmbedLength][sourceAgent]];
				sourceVal[c] *= base;
				sourceVal[c] += states[r][sourceAgent];
			}
		}		
	}

	/**
 	 * Add observations in to our estimates of the PDFs,
 	 * from a multivariate time-series.
 	 * This call suitable only for homogeneous agents, as all
 	 *  variable pairs separated by h rows and j columns
 	 *  will contribute to the PDFs.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param h - number of rows to compute transfer entropy across
	 * 	(i.e. source is in row i-h, dest is column i) 
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. source is column i-j, dest is column i) 
	 */
	public void addObservations(int states[][][], int h, int j) {
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
		if (periodicBoundaryConditions) {
			observations += (timeSteps - startObservationTime) * agentRows * agentColumns;
		} else {
			observations += (timeSteps - startObservationTime) * (agentRows - Math.abs(h)) * (agentColumns - Math.abs(j));
		}
		
		// Initialise and store the current previous and source value for each agent
		int[][] pastVal = new int[agentRows][agentColumns];
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
				pastVal[r][c] = 0;
				for (int p = 0; p < k; p++) {
					pastVal[r][c] *= base;
					pastVal[r][c] += states[startObservationTime - k + p][r][c];
				}
			}
		}
		int[][] sourceVal = new int[agentRows][agentColumns];
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
				sourceVal[r][c] = 0;
				int sourceAgentRow = r-h;
				if ((sourceAgentRow < 0) || (sourceAgentRow >= agentRows)) {
					// Source agent is out of bounds unless we are using periodic boundary conditions
					if (periodicBoundaryConditions) {
						sourceAgentRow = (sourceAgentRow+agentRows) % agentRows;
					} else {
						// Don't add this to our observations
						continue;
					}
				}
				int sourceAgentColumn = c-j;
				if ((sourceAgentColumn < 0) || (sourceAgentColumn >= agentColumns)) {
					// Source agent is out of bounds unless we are using periodic boundary conditions
					if (periodicBoundaryConditions) {
						sourceAgentColumn = (sourceAgentColumn+agentColumns) % agentColumns;
					} else {
						// Don't add this to our observations
						continue;
					}
				}
				for (int p = 0; p < sourceHistoryEmbedLength; p++) {
					sourceVal[r][c] *= base;
					sourceVal[r][c] += states[startObservationTime - sourceHistoryEmbedLength + p][sourceAgentRow][sourceAgentColumn];
				}
			}
		}

		// 1. Count the tuples observed
		int destVal;
		for (int t = startObservationTime; t < timeSteps; t++) {
			for (int r = 0; r < agentRows; r++) {
				for (int c = 0; c < agentColumns; c++) {
					// Add to the count for this particular transition:
					// (cell's assigned as above)
					int sourceAgentRow = r-h;
					if ((sourceAgentRow < 0) || (sourceAgentRow >= agentRows)) {
						// Source agent is out of bounds unless we are using periodic boundary conditions
						if (periodicBoundaryConditions) {
							sourceAgentRow = (sourceAgentRow+agentRows) % agentRows;
						} else {
							// Don't add this to our observations
							continue;
						}
					}
					int sourceAgentColumn = c-j;
					if ((sourceAgentColumn < 0) || (sourceAgentColumn >= agentColumns)) {
						// Source agent is out of bounds unless we are using periodic boundary conditions
						if (periodicBoundaryConditions) {
							sourceAgentColumn = (sourceAgentColumn+agentColumns) % agentColumns;
						} else {
							// Don't add this to our observations
							continue;
						}
					}
					destVal = states[t][r][c];
					sourceNextPastCount[sourceVal[r][c]][destVal][pastVal[r][c]]++;
					sourcePastCount[sourceVal[r][c]][pastVal[r][c]]++;
					nextPastCount[destVal][pastVal[r][c]]++;
					pastCount[pastVal[r][c]]++;
					nextCount[destVal]++;
					// Update the previous value:
					if (k > 0) {
						pastVal[r][c] -= maxShiftedValue[states[t-k][r][c]];
						pastVal[r][c] *= base;
						pastVal[r][c] += states[t][r][c];
					}
					sourceVal[r][c] -= maxShiftedSourceValue[states[t-sourceHistoryEmbedLength][sourceAgentRow][sourceAgentColumn]];
					sourceVal[r][c] *= base;
					sourceVal[r][c] += states[t][sourceAgentRow][sourceAgentColumn];
				}
			}
		}		
	}

	/**
 	 * Add observations for a single source-destination pair of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to {@link #addObservations(int[][], int)}
 	 *  for computing TE for heterogeneous agents.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param sourceIndex source variable index in states
	 * @param destIndex destination variable index in states
	 */
	public void addObservations(int states[][], int sourceIndex, int destIndex) {
		int rows = states.length;
		// increment the count of observations:
		observations += (rows - startObservationTime); 
		
		// Initialise and store the current previous value for each column
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += states[startObservationTime - k + p][destIndex];
		}
		int sourceVal = 0;
		for (int p = 0; p < sourceHistoryEmbedLength; p++) {
			sourceVal *= base;
			sourceVal += states[startObservationTime - sourceHistoryEmbedLength + p][sourceIndex];
		}

		// 1. Count the tuples observed
		int destVal;
		for (int r = startObservationTime; r < rows; r++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			destVal = states[r][destIndex];
			sourceNextPastCount[sourceVal][destVal][pastVal]++;
			sourcePastCount[sourceVal][pastVal]++;
			nextPastCount[destVal][pastVal]++;
			pastCount[pastVal]++;
			nextCount[destVal]++;
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[states[r-k][destIndex]];
				pastVal *= base;
				pastVal += states[r][destIndex];
			}
			sourceVal -= maxShiftedSourceValue[states[r-sourceHistoryEmbedLength][sourceIndex]];
			sourceVal *= base;
			sourceVal += states[r][sourceIndex];
		}
	}

	/**
 	 * Add observations for a single source-destination pair of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to {@link #addObservations(int[][][], int, int)}
 	 *  for computing TE for heterogeneous agents.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param sourceRowIndex source variable row index in states
	 * @param sourceColumnIndex source variable column index in states
	 * @param destRowIndex destination variable row index in states
	 * @param destColumnIndex destination variable column index in states
	 */
	public void addObservations(int states[][][], int sourceRowIndex, int sourceColumnIndex,
												  int destRowIndex, int destColumnIndex) {
		int timeSteps = states.length;
		// increment the count of observations:
		observations += (timeSteps - startObservationTime); 
		
		// Initialise and store the current previous value for each column
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += states[startObservationTime - k + p][destRowIndex][destColumnIndex];
		}
		int sourceVal = 0;
		for (int p = 0; p < sourceHistoryEmbedLength; p++) {
			sourceVal *= base;
			sourceVal += states[startObservationTime - sourceHistoryEmbedLength + p][sourceRowIndex][sourceColumnIndex];
		}
		
		// 1. Count the tuples observed
		int destVal;
		for (int r = startObservationTime; r < timeSteps; r++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			destVal = states[r][destRowIndex][destColumnIndex];
			sourceNextPastCount[sourceVal][destVal][pastVal]++;
			sourcePastCount[sourceVal][pastVal]++;
			nextPastCount[destVal][pastVal]++;
			pastCount[pastVal]++;
			nextCount[destVal]++;
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[states[r-k][destRowIndex][destColumnIndex]];
				pastVal *= base;
				pastVal += states[r][destRowIndex][destColumnIndex];
			}
			sourceVal -= maxShiftedSourceValue[states[r-sourceHistoryEmbedLength][sourceRowIndex][sourceColumnIndex]];
			sourceVal *= base;
			sourceVal += states[r][sourceRowIndex][sourceColumnIndex];
		}
	}

	/**
	 * 
	 * Returns the count of observations of the supplied past state
	 *  pastVal.
	 * The past state is indicated by a unique discrete integer representing the joint variable
	 *  of the k past states: (dest[n-k+1],dest[n-k+2],...,dest[n-1],dest[n]).
	 * The integer is computed as:<br/>
	 * pastVal = dest[n-k+1] * base^(k-1) + dest[n-k+2] * base^(k-2) + ... + dest[n-1] * base + dest[n]
	 * 
	 * 
	 * @param pastVal int representing the joint state of the past of the destination dest[n]^k
	 * @return count of observations of this given past state
	 */
	public int getPastCount(int pastVal) {
		return pastCount[pastVal];
	}
	
	/**
	 * Returns the probability of the supplied past state
	 *  pastVal.
	 * See {@link #getPastCount(int)} for how the joint value representing the past is calculated.
	 * 
	 * @param pastVal int representing the joint state of the past of the destination dest[n]^k
	 * @return probability of the given past state
	 */
	public double getPastProbability(int pastVal) {
		return (double) pastCount[pastVal] / (double) observations;
	}

	/**
	 * Returns the count of observations of the past given past state and next value.
	 * 
	 * See {@link #getPastCount(int)} for how the joint value representing the past is calculated.
	 * 
	 * @param destVal next state of the destination dest[n+1]
	 * @param pastVal int representing the joint state of the past of the destination dest[n]^k
	 * @return count of observations of the given past state and next state
	 */
	public int getNextPastCount(int destVal, int pastVal) {
		return nextPastCount[destVal][pastVal];
	}
	
	/**
	 * Returns the probability of the past given past state and next value.
	 * 
	 * See {@link #getPastCount(int)} for how the joint value representing the past is calculated.
	 * 
	 * @param destVal next state of the destination dest[n+1]
	 * @param pastVal int representing the joint state of the past of the destination dest[n]^k
	 * @return probability of the given past state and next state
	 */
	public double getNextPastProbability(int destVal, int pastVal) {
		return (double) nextPastCount[destVal][pastVal] / (double) observations;
	}

	/**
	 * Returns the count of observations of the past given state dest[n]^k and the source state source[n]^l.
	 * 
	 * See {@link #getPastCount(int)} for how the joint values representing the past states are calculated.
	 * 
	 * @param sourceVal int representing the joint state of the source source[n]^l
	 * @param pastVal int representing the joint state of the past of the destination dest[n]^k
	 * @return count of observations of the given past state and the source state
	 */
	public int getSourcePastCount(int sourceVal, int pastVal) {
		return sourcePastCount[sourceVal][pastVal];
	}
	
	/**
	 * Returns the probability of the past given state dest[n]^k and the source state source[n]^l.
	 * 
	 * See {@link #getPastCount(int)} for how the joint values representing the past states are calculated.
	 * 
	 * @param sourceVal int representing the joint state of the source source[n]^l
	 * @param pastVal int representing the joint state of the past of the destination dest[n]^k
	 * @return probability of the given past state and the source state
	 */
	public double getSourcePastProbability(int sourceVal, int pastVal) {
		return (double) sourcePastCount[sourceVal][pastVal] / (double) observations;
	}

	/**
	 * Returns the count of observations of the past given state dest[n]^k,
	 *  the next state of the destination dest[n+1] and the source state source[n]^l.
	 *  
	 * See {@link #getPastCount(int)} for how the joint values representing the past states are calculated.
	 * 
	 * @param sourceVal int representing the joint state of the source source[n]^l
	 * @param nextVal next state of the destination dest[n+1]
	 * @param pastVal int representing the joint state of the past of the destination dest[n]^k
	 * @return count of observations of the given past state, next state of destination and the source state
	 */
	public int getSourceNextPastCount(int sourceVal, int destVal, int pastVal) {
		return sourceNextPastCount[sourceVal][destVal][pastVal];
	}
	
	/**
	 * Returns the probability of the past given state dest[n]^k,
	 *  the next state of the destination dest[n+1] and the source state source[n]^l.
	 *  
	 * See {@link #getPastCount(int)} for how the joint values representing the past states are calculated.
	 * 
	 * @param sourceVal int representing the joint state of the source source[n]^l
	 * @param nextVal next state of the destination dest[n+1]
	 * @param pastVal int representing the joint state of the past of the destination dest[n]^k
	 * @return probability of the given past state, next state of destination and the source state
	 */
	public double getSourceNextPastProbability(int sourceVal, int destVal, int pastVal) {
		return (double) sourceNextPastCount[sourceVal][destVal][pastVal] / (double) observations;
	}

	/**
	 * Returns the count of observations of the next state dest[n+1].
	 * 
	 * @param nextVal next state of the destination dest[n+1]
	 * @return count of observations of the given next state
	 */
	public int getNextCount(int destVal) {
		return nextCount[destVal];
	}
	
	/**
	 * Returns the probability of the next state dest[n+1].
	 * 
	 * @param nextVal state of the next destination dest[n+1]
	 * @return probability of the given next state
	 */
	public double getNextProbability(int destVal) {
		return (double) nextCount[destVal] / (double) observations;
	}

	@Override
	public double computeAverageLocalOfObservations() {
		double te = 0.0;
		double teCont = 0.0;

		max = 0;
		min = 0;
		double meanSqLocals = 0;
		for (int pastVal = 0; pastVal < base_power_k; pastVal++) {
			// compute p(past)
			// double p_past = (double) pastCount[pastVal] / (double) observations;
			for (int destVal = 0; destVal < base; destVal++) {
				// compute p(dest,past)
				// double p_dest_past = (double) destPastCount[destVal][pastVal] / (double) observations;
				for (int sourceVal = 0; sourceVal < base_power_l; sourceVal++) {
					// compute p(source,dest,past)
					double p_source_dest_past = (double) sourceNextPastCount[sourceVal][destVal][pastVal] / (double) observations;
					// compute p(source,past)
					// double p_source_past = (double) sourcePastCount[sourceVal][pastVal] / (double) observations;
					// Compute TE contribution:
					if (sourceNextPastCount[sourceVal][destVal][pastVal] != 0) {
						/* Double check: should never happen
						if ((sourcePastCount[sourceVal][pastVal] == 0) ||
							(destPastCount[destVal][pastVal] == 0) ||
							(pastCount[pastVal] == 0)) {
							throw new RuntimeException("one subcount was zero!!");
						}
						*/
						
						double logTerm = ((double) sourceNextPastCount[sourceVal][destVal][pastVal] / (double) sourcePastCount[sourceVal][pastVal]) /
						 	((double) nextPastCount[destVal][pastVal] / (double) pastCount[pastVal]);
						double localValue = Math.log(logTerm) / log_2;
						teCont = p_source_dest_past * localValue;
						if (localValue > max) {
							max = localValue;
						} else if (localValue < min) {
							min = localValue;
						}
						// Add this contribution to the mean 
						//  of the squared local values
						meanSqLocals += teCont * localValue;
					} else {
						teCont = 0.0;
					}
					te += teCont;
				}
			}
		}
		
		average = te;
		std = Math.sqrt(meanSqLocals - average * average);
		estimateComputed = true;
		return te;
	}
	
	/**
	 * Returns the average active information storage from
	 *  the observed values which have been passed in previously. 
	 *  
	 * @see ActiveInformationCalculatorDiscrete
	 */
	public double computeAverageActiveInfoStorageOfObservations() {
		double active = 0.0;
		double activeCont = 0.0;

		for (int nextVal = 0; nextVal < base; nextVal++) {
			// compute p_next
			double p_next = (double) nextCount[nextVal] / (double) observations;
			for (int prevVal = 0; prevVal < base_power_k; prevVal++) {
				// Compute MI contribution:
				if (nextPastCount[nextVal][prevVal] != 0) {
					double logTerm = (double) nextPastCount[nextVal][prevVal] /
								(double) pastCount[prevVal] /
								p_next;
					double localValue = Math.log(logTerm) / log_2;
					activeCont = (nextPastCount[nextVal][prevVal] /
								(double) observations) * localValue;
				} else {
					activeCont = 0.0;
				}
				active += activeCont;
			}
		}
		
		return active;
	}

	/**
	 * Dump a debug print of the PDFs of our observations
	 */
	public void debugPrintObservations() {
		
		System.out.println("Src\tDst\tPast\tc(s,d,p)\tc(s,p)\tc(d,p)\tc(p)");
		for (int pastVal = 0; pastVal < base_power_k; pastVal++) {
			for (int destVal = 0; destVal < base; destVal++) {
				for (int sourceVal = 0; sourceVal < base_power_l; sourceVal++) {
					// Compute TE contribution:
					System.out.println(sourceVal + "\t" + destVal + "\t" + pastVal + "\t" +
							sourceNextPastCount[sourceVal][destVal][pastVal] + "\t\t" +
							sourcePastCount[sourceVal][pastVal] + "\t" + 
							nextPastCount[destVal][pastVal] + "\t" +
							pastCount[pastVal]);
				}
			}
		}
	}

	@Override
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) {
		double actualTE = computeAverageLocalOfObservations();
		
		// Reconstruct the *joint* source values (not necessarily in order, but using joint values retains their l-tuples)
		int[] sourceValues = new int[observations];
		int t_s = 0;
		for (int sourceVal = 0; sourceVal < base_power_l; sourceVal++) {
			// Count up the number of times this joint source value was observed:
			int numberOfSamples = 0;
			for (int pastVal = 0; pastVal < base_power_k; pastVal++) {
				numberOfSamples += sourcePastCount[sourceVal][pastVal];
			}
			// Now add all of these as unordered observations:
			MatrixUtils.fill(sourceValues, sourceVal, t_s, numberOfSamples);
			t_s += numberOfSamples;
		}
		
		// And construct unordered (dest,past) tuples.
		// It doesn't matter that we've appeared to destroy the ordering here because
		//  the joint distribution nextPastCount is actually preserved in
		//  our construction of pastVal and destValues together here.
		int[] destValues = new int[observations];
		int[] pastValues = new int[observations];
		int t_d = 0;
		int t_p = 0;
		for (int pastVal = 0; pastVal < base_power_k; pastVal++) {
			MatrixUtils.fill(pastValues, pastVal, t_p, pastCount[pastVal]);
			t_p += pastCount[pastVal];
			for (int destVal = 0; destVal < base; destVal++) {
				MatrixUtils.fill(destValues, destVal, t_d, nextPastCount[destVal][pastVal]);
				t_d += nextPastCount[destVal][pastVal];
			}
		}
		
		// Construct new source orderings based on the source probabilities only
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(observations, numPermutationsToCheck);

		// TODO The use of base_power_l as the base for all variables is particularly wasteful
		//  of resources here, but there's not much other choice. A better solution
		//  will come when we switch to an underlying conditional MI calculator, with separate
		//  bases for each variable, and just use its computeSignificance() method.
		TransferEntropyCalculatorDiscrete ate2 = new TransferEntropyCalculatorDiscrete(base_power_l, k, 1);
		ate2.initialise();
		ate2.observations = observations;
		ate2.pastCount = pastCount;
		ate2.nextPastCount = nextPastCount;
		int countWhereTeIsMoreSignificantThanOriginal = 0;
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		for (int p = 0; p < numPermutationsToCheck; p++) {
			// Generate a new re-ordered data set for the source
			int[] newSourceData = MatrixUtils.extractSelectedTimePoints(sourceValues, newOrderings[p]);
			// compute the joint probability distributions
			MatrixUtils.fill(ate2.sourceNextPastCount, 0);
			MatrixUtils.fill(ate2.sourcePastCount, 0);
			for (int t = 0; t < observations; t++) {
				// This looks like we're doing no source embedding --
				//  but we're already using joint source values in newSourceData
				ate2.sourcePastCount[newSourceData[t]][pastValues[t]]++;
				ate2.sourceNextPastCount[newSourceData[t]][destValues[t]][pastValues[t]]++;
			}
			// And get a TE value for this realisation:
			double newTe = ate2.computeAverageLocalOfObservations();
			measDistribution.distribution[p] = newTe;
			if (newTe >= actualTE) {
				countWhereTeIsMoreSignificantThanOriginal++;
			}

		}
		
		// And return the significance
		measDistribution.pValue = (double) countWhereTeIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualTE;
		return measDistribution;
	}
	
	@Override
	public AnalyticMeasurementDistribution computeSignificance()
			throws Exception {
		if (!estimateComputed) {
			computeAverageLocalOfObservations();
		}
		return new ChiSquareMeasurementDistribution(2.0*((double)observations)*average,
				(sourceHistoryEmbedLength*base - 1)*(base - 1)*(k*base));
	}

	/**
	 * Computes local transfer entropy for the given values
	 * 
	 * See {@link #getPastCount(int)} for how the joint values representing the past are calculated.
	 * 
	 * @param sourceCurrent int representing the joint state of the source source[n]^l
	 * @param destNext next state of the destination dest[n+1]
	 * @param destPast int representing the joint state of the past of the destination dest[n]^k
	 * 
	 * @return local TE for the given observation
	 */
	public double computeLocalFromPreviousObservations(int sourceCurrent, int destNext, int destPast){

		double logTerm = ((double) sourceNextPastCount[sourceCurrent][destNext][destPast] / (double) sourcePastCount[sourceCurrent][destPast]) /
	 		((double) nextPastCount[destNext][destPast] / (double) pastCount[destPast]);
		return Math.log(logTerm) / log_2;
	}

	/**
	 * Computes local apparent transfer entropy for the given
	 *  states, using PDFs built up from observations previously
	 *  sent in via the addObservations method.
	 *  
 	 * @param source source time-series
	 * @param dest destination time-series. 
	 *  Must be same length as source
	 * @return time-series of local TE values
	 */
	public double[] computeLocalFromPreviousObservations(int sourceStates[], int destStates[]){
		int timeSteps = destStates.length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localTE = new double[timeSteps];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int pastVal = 0; 
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += destStates[startObservationTime - k + p];
		}
		int sourceVal = 0;
		for (int p = 0; p < sourceHistoryEmbedLength; p++) {
			sourceVal *= base;
			sourceVal += sourceStates[startObservationTime - sourceHistoryEmbedLength + p];
		}
		int destVal;
		double logTerm;
		for (int t = startObservationTime; t < timeSteps; t++) {
			destVal = destStates[t];
			// Now compute the local value
			logTerm = ((double) sourceNextPastCount[sourceVal][destVal][pastVal] / (double) sourcePastCount[sourceVal][pastVal]) /
		 		((double) nextPastCount[destVal][pastVal] / (double) pastCount[pastVal]);
			localTE[t] = Math.log(logTerm) / log_2;
			average += localTE[t];
			if (localTE[t] > max) {
				max = localTE[t];
			} else if (localTE[t] < min) {
				min = localTE[t];
			}
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[destStates[t-k]];
				pastVal *= base;
				pastVal += destStates[t];
			}
			sourceVal -= maxShiftedSourceValue[sourceStates[t-sourceHistoryEmbedLength]];
			sourceVal *= base;
			sourceVal += sourceStates[t];
		}		

		average = average/(double) (timeSteps - startObservationTime);
		
		return localTE;
	}

	/**
	 * Computes local transfer for the given
	 *  multivariate states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
 	 * This call suitable only for homogeneous agents, as all
 	 *  variable pairs separated by j column will
 	 *  have their local TE computed.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. source is column i-j, dest is column i: we
	 *  compute transfer is j cells to the right, using observations
	 *  across all column pairs separated by j) 
	 * @return multivariate time series of local TE values
	 *  (first index is time, second index is destination variable)
	 */
	public double[][] computeLocalFromPreviousObservations(int states[][], int j){
		int timeSteps = states.length;
		int agents = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][] localTE = new double[timeSteps][agents];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int[] pastVal = new int[agents]; 
		for (int c = 0; c < agents; c++) {
			pastVal[c] = 0;
			for (int p = 0; p < k; p++) {
				pastVal[c] *= base;
				pastVal[c] += states[startObservationTime - k + p][c];
			}
		}
		int[] sourceVal = new int[agents];
		for (int c = 0; c < agents; c++) {
			sourceVal[c] = 0;
			int sourceAgentIndex = c-j;
			if ((sourceAgentIndex < 0) || (sourceAgentIndex >= agents)) {
				// Source agent is out of bounds unless we are using periodic boundary conditions
				if (periodicBoundaryConditions) {
					sourceAgentIndex = (sourceAgentIndex+agents) % agents;
				} else {
					// Don't add this to our observations
					continue;
				}
			}
			for (int p = 0; p < sourceHistoryEmbedLength; p++) {
				sourceVal[c] *= base;
				sourceVal[c] += states[startObservationTime - sourceHistoryEmbedLength + p][sourceAgentIndex];
			}
		}
		int destVal;
		double logTerm;
		for (int t = startObservationTime; t < timeSteps; t++) {
			for (int c = 0; c < agents; c++) {
				int sourceAgentIndex = c-j;
				if ((sourceAgentIndex < 0) || (sourceAgentIndex >= agents)) {
					// Source agent is out of bounds unless we are using periodic boundary conditions
					if (periodicBoundaryConditions) {
						sourceAgentIndex = (sourceAgentIndex+agents) % agents;
					} else {
						// Don't compute a local value for this one
						continue;
					}
				}
				destVal = states[t][c];
				// Now compute the local value
				logTerm = ((double) sourceNextPastCount[sourceVal[c]][destVal][pastVal[c]] / (double) sourcePastCount[sourceVal[c]][pastVal[c]]) /
			 		((double) nextPastCount[destVal][pastVal[c]] / (double) pastCount[pastVal[c]]);
				localTE[t][c] = Math.log(logTerm) / log_2;
				average += localTE[t][c];
				if (localTE[t][c] > max) {
					max = localTE[t][c];
				} else if (localTE[t][c] < min) {
					min = localTE[t][c];
				}
				// Update the previous value:
				if (k > 0) {
					pastVal[c] -= maxShiftedValue[states[t-k][c]];
					pastVal[c] *= base;
					pastVal[c] += states[t][c];
				}
				sourceVal[c] -= maxShiftedSourceValue[states[t-sourceHistoryEmbedLength][sourceAgentIndex]];
				sourceVal[c] *= base;
				sourceVal[c] += states[t][sourceAgentIndex];
			}
		}		

		if (periodicBoundaryConditions) {
			average = average/(double) ((timeSteps - startObservationTime) * agents);
		} else {
			average = average/(double) ((timeSteps - startObservationTime) * (agents - Math.abs(j)));
		}
		
		return localTE;
	}
	
	/**
	 * Computes local transfer for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
 	 * This call suitable only for homogeneous agents, as all
 	 *  variable pairs separated by h rows and j columns
 	 *  will have their local TE computed.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param h - number of rows to compute transfer entropy across
	 * 	(i.e. source is in row i-h, dest is column i) 
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. source is column i-j, dest is column i) 
	 * @return multivariate time series of local TE values
	 *  (first index is time, second index is destination variable
	 *  row number, third is destination variable column number)
	 */
	public double[][][] computeLocalFromPreviousObservations(int states[][][], int h, int j){
		int timeSteps = states.length;
		int agentRows = states[0].length;
		int agentColumns = states[0][0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][][] localTE = new double[timeSteps][agentRows][agentColumns];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int[][] pastVal = new int[agentRows][agentColumns]; 
		for (int r = 0; r < agentRows; r++){
			for (int c = 0; c < agentColumns; c++) {
				pastVal[r][c] = 0;
				for (int p = 0; p < k; p++) {
					pastVal[r][c] *= base;
					pastVal[r][c] += states[startObservationTime - k + p][r][c];
				}
			}
		}
		int[][] sourceVal = new int[agentRows][agentColumns];
		for (int r = 0; r < agentRows; r++) {
			for (int c = 0; c < agentColumns; c++) {
				sourceVal[r][c] = 0;
				int sourceAgentRow = r-h;
				if ((sourceAgentRow < 0) || (sourceAgentRow >= agentRows)) {
					// Source agent is out of bounds unless we are using periodic boundary conditions
					if (periodicBoundaryConditions) {
						sourceAgentRow = (sourceAgentRow+agentRows) % agentRows;
					} else {
						// Don't add this to our observations
						continue;
					}
				}
				int sourceAgentColumn = c-j;
				if ((sourceAgentColumn < 0) || (sourceAgentColumn >= agentColumns)) {
					// Source agent is out of bounds unless we are using periodic boundary conditions
					if (periodicBoundaryConditions) {
						sourceAgentColumn = (sourceAgentColumn+agentColumns) % agentColumns;
					} else {
						// Don't add this to our observations
						continue;
					}
				}
				for (int p = 0; p < sourceHistoryEmbedLength; p++) {
					sourceVal[r][c] *= base;
					sourceVal[r][c] += states[startObservationTime - sourceHistoryEmbedLength + p][sourceAgentRow][sourceAgentColumn];
				}
			}
		}

		int destVal;
		double logTerm;
		for (int t = startObservationTime; t < timeSteps; t++) {
			for (int r = 0; r < agentRows; r++) {
				for (int c = 0; c < agentColumns; c++) {
					int sourceAgentRow = r-h;
					if ((sourceAgentRow < 0) || (sourceAgentRow >= agentRows)) {
						// Source agent is out of bounds unless we are using periodic boundary conditions
						if (periodicBoundaryConditions) {
							sourceAgentRow = (sourceAgentRow+agentRows) % agentRows;
						} else {
							// Don't compute a local value for this one
							continue;
						}
					}
					int sourceAgentColumn = c-j;
					if ((sourceAgentColumn < 0) || (sourceAgentColumn >= agentColumns)) {
						// Source agent is out of bounds unless we are using periodic boundary conditions
						if (periodicBoundaryConditions) {
							sourceAgentColumn = (sourceAgentColumn+agentColumns) % agentColumns;
						} else {
							// Don't compute a local value for this one
							continue;
						}
					}
					destVal = states[t][r][c];
					// Now compute the local value
					logTerm = ((double) sourceNextPastCount[sourceVal[r][c]][destVal][pastVal[r][c]] / (double) sourcePastCount[sourceVal[r][c]][pastVal[r][c]]) /
				 		((double) nextPastCount[destVal][pastVal[r][c]] / (double) pastCount[pastVal[r][c]]);
					localTE[t][r][c] = Math.log(logTerm) / log_2;
					average += localTE[t][r][c];
					if (localTE[t][r][c] > max) {
						max = localTE[t][r][c];
					} else if (localTE[t][r][c] < min) {
						min = localTE[t][r][c];
					}
					// Update the previous value:
					if (k > 0) {
						pastVal[r][c] -= maxShiftedValue[states[t-k][r][c]];
						pastVal[r][c] *= base;
						pastVal[r][c] += states[t][r][c];
					}
					sourceVal[r][c] -= maxShiftedSourceValue[states[t-sourceHistoryEmbedLength][sourceAgentRow][sourceAgentColumn]];
					sourceVal[r][c] *= base;
					sourceVal[r][c] += states[t][sourceAgentRow][sourceAgentColumn];
				}
			}
		}		

		if (periodicBoundaryConditions) {
			average = average/(double) ((timeSteps - startObservationTime) * agentRows * agentColumns);
		} else {
			average = average/(double) ((timeSteps - startObservationTime) * (agentRows - Math.abs(h)) *
												(agentColumns - Math.abs(j)));			
		}
		
		return localTE;
	}

	/**
	 * Computes local transfer for the given
	 *  single source-destination pair of the multi-agent system,
	 *  using pdfs built up from observations previously
	 *  sent in via the addObservations methods.
 	 * This call should be made as opposed to {@link #addObservations(int[][], int)}
 	 *  for computing local TE for heterogeneous agents.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param sourceIndex source variable index in states
	 * @param destIndex destination variable index in states
	 * @return time-series of local TE values between the series
	 */
	public double[] computeLocalFromPreviousObservations(int states[][], int sourceCol, int destCol){
		int rows = states.length;
		// int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localTE = new double[rows];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int pastVal = 0; 
		pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += states[startObservationTime - k + p][destCol];
		}
		int sourceVal = 0;
		for (int p = 0; p < sourceHistoryEmbedLength; p++) {
			sourceVal *= base;
			sourceVal += states[startObservationTime - sourceHistoryEmbedLength + p][sourceCol];
		}
		int destVal;
		double logTerm;
		for (int r = startObservationTime; r < rows; r++) {
			destVal = states[r][destCol];
			// Now compute the local value
			logTerm = ((double) sourceNextPastCount[sourceVal][destVal][pastVal] / (double) sourcePastCount[sourceVal][pastVal]) /
		 		((double) nextPastCount[destVal][pastVal] / (double) pastCount[pastVal]);
			localTE[r] = Math.log(logTerm) / log_2;
			average += localTE[r];
			if (localTE[r] > max) {
				max = localTE[r];
			} else if (localTE[r] < min) {
				min = localTE[r];
			}
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[states[r-k][destCol]];
				pastVal *= base;
				pastVal += states[r][destCol];
			}
			sourceVal -= maxShiftedSourceValue[states[r-sourceHistoryEmbedLength][sourceCol]];
			sourceVal *= base;
			sourceVal += states[r][sourceCol];
		}

		average = average/(double) (rows - startObservationTime);
		
		return localTE;
	}

	/**
	 * Computes local transfer for the given
	 *  single source-destination pair of the multi-agent system,
	 *  using pdfs built up from observations previously
	 *  sent in via the addObservations method.
 	 * This call should be made as opposed to {@link #addObservations(int[][][], int, int)}
 	 *  for computing local TE for heterogeneous agents.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param sourceRowIndex source variable row index in states
	 * @param sourceColumnIndex source variable column index in states
	 * @param destRowIndex destination variable row index in states
	 * @param destColumnIndex destination variable column index in states
	 * @return time-series of local TE values between the series
	 */
	public double[] computeLocalFromPreviousObservations(int states[][][],
			int sourceRowIndex, int sourceColumnIndex, int destRowIndex, int destColumnIndex){
		int timeSteps = states.length;
		// int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localTE = new double[timeSteps];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int pastVal = 0; 
		pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += states[startObservationTime - k + p][destRowIndex][destColumnIndex];
		}
		int sourceVal = 0;
		for (int p = 0; p < sourceHistoryEmbedLength; p++) {
			sourceVal *= base;
			sourceVal += states[startObservationTime - sourceHistoryEmbedLength + p][sourceRowIndex][sourceColumnIndex];
		}
		int destVal;
		double logTerm;
		for (int r = startObservationTime; r < timeSteps; r++) {
			destVal = states[r][destRowIndex][destColumnIndex];
			// Now compute the local value
			logTerm = ((double) sourceNextPastCount[sourceVal][destVal][pastVal] / (double) sourcePastCount[sourceVal][pastVal]) /
		 		((double) nextPastCount[destVal][pastVal] / (double) pastCount[pastVal]);
			localTE[r] = Math.log(logTerm) / log_2;
			average += localTE[r];
			if (localTE[r] > max) {
				max = localTE[r];
			} else if (localTE[r] < min) {
				min = localTE[r];
			}
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[states[r-k][destRowIndex][destColumnIndex]];
				pastVal *= base;
				pastVal += states[r][destRowIndex][destColumnIndex];
			}
			sourceVal -= maxShiftedSourceValue[states[r-sourceHistoryEmbedLength][sourceRowIndex][sourceColumnIndex]];
			sourceVal *= base;
			sourceVal += states[r][sourceRowIndex][sourceColumnIndex];
		}

		average = average/(double) (timeSteps - startObservationTime);
		
		return localTE;
	}

	/**
	 * Standalone routine to 
	 * compute local transfer entropy between two time series
	 * Return a time series of local values.
	 * First max(k,l) values are zeros since TE is not defined there
	 * 
 	 * @param sourceStates source time-series
	 * @param destStates destination time-series. 
	 *  Must be same length as sourceStates
	 * @return time-series of local TE values
	 */
	public double[] computeLocal(int sourceStates[], int destStates[]) {
		
		initialise();
		addObservations(sourceStates, destStates);
		return computeLocalFromPreviousObservations(sourceStates, destStates);
	}

	/**
	 * Standalone routine to 
	 * compute local transfer entropy across a 2D spatiotemporal
	 *  array of the states of homogeneous agents.
	 * Return a 2D spatiotemporal array of local values.
	 * First max(k,l) values are zeros since TE is not defined there.
 	 * This call suitable only for homogeneous agents, as all
 	 *  variable pairs separated by j column will
 	 *  have their local TE computed.
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. source is column i-j, dest is column i: we
	 *  compute transfer is j cells to the right, using observations
	 *  across all column pairs separated by j) 
	 * @return multivariate time series of local TE values
	 *  (first index is time, second index is destination variable)
	 */
	public double[][] computeLocal(int states[][], int j) {
		
		initialise();
		addObservations(states, j);
		return computeLocalFromPreviousObservations(states, j);
	}
	
	/**
	 * Standalone routine to 
	 * compute local transfer entropy across a 3D spatiotemporal
	 *  array of the states of homogeneous agents.
	 * Return a 3D spatiotemporal array of local values.
	 * First max(k,l) values are zeros since TE is not defined there.
 	 * This call suitable only for homogeneous agents, as all
 	 *  variable pairs separated by h rows and j columns
 	 *  will have their local TE computed.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param h - number of rows to compute transfer entropy across
	 * 	(i.e. source is in row i-h, dest is column i) 
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. source is column i-j, dest is column i) 
	 * @return multivariate time series of local TE values
	 *  (first index is time, second index is destination variable
	 *  row number, third is destination variable column number)
	 */
	public double[][][] computeLocal(int states[][][], int h, int j) {
		
		initialise();
		addObservations(states, h, j);
		return computeLocalFromPreviousObservations(states, h, j);
	}

	/**
	 * Standalone routine to 
	 * compute average local transfer entropy across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return the average TE.
 	 * This call suitable only for homogeneous agents, as all
 	 *  variable pairs separated by j column will
 	 *  have their local TE computed.
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. source is column i-j, dest is column i: we
	 *  compute transfer is j cells to the right, using observations
	 *  across all column pairs separated by j) 
	 * @return average TE across j variables to the right
	 */
	public double computeAverageLocal(int states[][], int j) {
		
		initialise();
		addObservations(states, j);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute average local transfer entropy across a 3D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return the average.
 	 * This call suitable only for homogeneous agents, as all
 	 *  variable pairs separated by h rows and j columns
 	 *  will have their PDFs combined.
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param h - number of rows to compute transfer entropy across
	 * 	(i.e. source is in row i-h, dest is column i) 
	 * @param j - number of columns to compute transfer entropy across
	 * 	(i.e. source is column i-j, dest is column i) 
	 * @return
	 */
	public double computeAverageLocal(int states[][][], int h, int j) {
		
		initialise();
		addObservations(states, h, j);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute local transfer entropy between specific variables in
	 * a 2D spatiotemporal multivariate time-series.
	 * First max(k,l) values are zeros since TE is not defined there.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param sourceCol source variable index in states
	 * @param destCol destination variable index in states
	 * @return time-series of local TE values between the series
	 */
	public double[] computeLocal(int states[][], int sourceCol, int destCol) {
		
		initialise();
		addObservations(states, sourceCol, destCol);
		return computeLocalFromPreviousObservations(states, sourceCol, destCol);
	}

	/**
	 * Standalone routine to 
	 * computes local transfer for the given
	 *  single source-destination pair of the 3D multi-agent system.
	 * This method suitable for heterogeneous variables.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param sourceRowIndex source variable row index in states
	 * @param sourceColumnIndex source variable column index in states
	 * @param destRowIndex destination variable row index in states
	 * @param destColumnIndex destination variable column index in states
	 * @return time-series of local TE values between the series
	 */
	public double[] computeLocal(int states[][][], int sourceRowIndex, int sourceColumnIndex,
											int destRowIndex, int destColumnIndex) {
		
		initialise();
		addObservations(states, sourceRowIndex, sourceColumnIndex, destRowIndex, destColumnIndex);
		return computeLocalFromPreviousObservations(states, sourceRowIndex, sourceColumnIndex,
				destRowIndex, destColumnIndex);
	}

	/**
	 * Standalone routine to 
	 * compute local transfer entropy between specific variables in
	 * a 2D spatiotemporal multivariate time-series.
	 * Returns the average.
	 * This method suitable for heterogeneous agents.
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param sourceCol source variable index in states
	 * @param destCol destination variable index in states
	 * @return average TE for the given pair
	 */
	public double computeAverageLocal(int states[][], int sourceCol, int destCol) {
		
		initialise();
		addObservations(states, sourceCol, destCol);
		return computeAverageLocalOfObservations();
	}
	
	/**
	 * Standalone routine to 
	 * compute local transfer entropy between specific variables in
	 * a 3D spatiotemporal multivariate time-series.
	 * Returns the average.
	 * This method suitable for heterogeneous agents.
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param sourceRowIndex source variable row index in states
	 * @param sourceColumnIndex source variable column index in states
	 * @param destRowIndex destination variable row index in states
	 * @param destColumnIndex destination variable column index in states
	 * @return average TE for the given pair
	 */
	public double computeAverageLocal(int states[][][], int sourceRowIndex, int sourceColumnIndex,
			int destRowIndex, int destColumnIndex) {
		
		initialise();
		addObservations(states, sourceRowIndex, sourceColumnIndex, destRowIndex, destColumnIndex);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Whether we assume periodic boundary conditions in the calls
	 *  for homogeneous variables.
	 *  
	 * @return as above
	 */
	public boolean isPeriodicBoundaryConditions() {
		return periodicBoundaryConditions;
	}
	/**
	 * set whether we assume periodic boundary conditions in the calls
	 *  for homogeneous variables.
	 *  
	 * @param periodicBoundaryConditions as above
	 */
	public void setPeriodicBoundaryConditions(boolean periodicBoundaryConditions) {
		this.periodicBoundaryConditions = periodicBoundaryConditions;
	}
}
