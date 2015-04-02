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
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

/**
 * <p>Implements <b>conditional transfer entropy</b>
 * for univariate discrete time-series data.
 * That is, it is applied to <code>int[]</code> data, indexed
 * by time.
 * See Schreiber below for the definition of transfer entropy,
 * and Lizier et al (2008, 2010). for the definition of local transfer entropy
 * and conditional TE, which is TE conditioned on one or 
 * more other potential sources.
 * This is also called complete TE when all other causal sources
 * are conditioned on.
 * </p>
 * 
 * <p>
 * Usage of the child classes implementing this interface is intended to follow this paradigm:
 * </p>
 * <ol>
 * 		<li>Construct the calculator via
 * 			{@link #ConditionalTransferEntropyCalculatorDiscrete(int, int, int)};</li>
 *		<li>Initialise the calculator using
 *			{@link #initialise()};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using one or more calls to
 * 			the set of {@link #addObservations(int[], int[], int[])} methods, then</li>
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
 * <p>
 * The conditional sources (specified using either their
 *  offsets from the destination variable or their absolute column numbers
 *  in the multivariate data set)
 *  should be supplied in the same order in every method call, otherwise the answer supplied will
 *  be incorrect.
 * </p>
 * 
 * <p><i>Note for developers</i>: Ideally, this class would extend ContextOfPastMeasure, however
 *  by conditioning on other info contributors, we need to alter
 *  the arrays pastCount and nextPastCount to consider all
 *  conditioned variables (i.e. other sources) also.
 * </p>
 * 
 * TODO Add methods for passing in single time series.
 * This is done for addObservations, but not other routines.
 * 
 * TODO Implement AnalyticNullDistributionComputer
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
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href=http://dx.doi.org/10.1063/1.3486801">
 *  "Information modification and particle collisions in distributed computation"</a>
 *  Chaos 20, 3, 037109 (2010).</li>
 * </ul>
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class ConditionalTransferEntropyCalculatorDiscrete extends InfoMeasureCalculatorDiscrete {

	protected int k = 0; // history length k.
	protected int base_power_k = 0;
	protected int base_power_num_others = 0;
	protected int numOtherInfoContributors = 0;
	protected int[][][][] sourceDestPastOthersCount = null;	// count for (i-j[n],i[n+1],i[n]^k,others) tuples
	protected int[][][] sourcePastOthersCount = null;			// count for (i-j[n],i[n]^k,others) tuples
	protected int[][][] destPastOthersCount = null; // Count for (i[n+1], i[n]^k,others) tuples
	protected int[][] pastOthersCount = null; // Count for (i[n]^k,others)
	protected int[] maxShiftedValue = null; // states * (base^(k-1))

	/**
	 * First time step at which we can take an observation
	 *  (needs to account for k previous steps)
	 */
	protected int startObservationTime = 1;
	
	/**
	 * User was formerly forced to create new instances through this factory method.
	 * Retained for backwards compatibility.
	 * 
	 * @param base
	 * @param history
	 * @param numOtherInfoContributors
	 * @deprecated
	 * @return
	 */
	public static ConditionalTransferEntropyCalculatorDiscrete
		newInstance(int base, int history, int numOtherInfoContributors) {
		
		return new ConditionalTransferEntropyCalculatorDiscrete
					(base, history, numOtherInfoContributors);
		
		// Old code for an attempted optimisation:
		/*
		if (isPowerOf2(base)) {
			return new CompleteTransferEntropyCalculatorBase2
						(base, history, numOtherInfoContributors);
		} else {
			return new CompleteTransferEntropyCalculator
					(base, history, numOtherInfoContributors);
		}
		*/
	}
	
	/**
	 * Construct a new instance
	 * 
	 * @param base number of symbols for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param history embedded history length of the destination to condition on -
	 *        this is k in Schreiber's notation.
	 * @param numOtherInfoContributors number of information contributors
	 *   (other than the past of the destination
	 *   or the source) to condition on.
	 */
	public ConditionalTransferEntropyCalculatorDiscrete
		(int base, int history, int numOtherInfoContributors) {

		super(base);
		
		k = history;
		this.numOtherInfoContributors = numOtherInfoContributors;
		base_power_k = MathsUtils.power(base, k);
		base_power_num_others = MathsUtils.power(base, numOtherInfoContributors);
		
		// Relaxing this assumption so we can use this calculation as
		//  a time-lagged conditional MI at will:
		//if (k < 1) {
		//	throw new RuntimeException("History k " + history + " is not >= 1 a ContextOfPastMeasureCalculator");
		//}
		
		// Which time step do we start taking observations from?
		// Normally this is k (to allow k previous time steps)
		//  but if k==0 (becoming a lagged MI), it's 1.
		startObservationTime = Math.max(k, 1);

		// check that we can convert the base tuple into an integer ok
		if (k > Math.log(Integer.MAX_VALUE) / log_base) {
			throw new RuntimeException("Base and history combination too large");
		}
		if (numOtherInfoContributors < 1) {
			throw new RuntimeException("Number of other info contributors < 1 for CompleteTECalculator");
		}
		
		// Create storage for counts of observations
		sourceDestPastOthersCount = new int[base][base][base_power_k][base_power_num_others];
		sourcePastOthersCount = new int[base][base_power_k][base_power_num_others];
		destPastOthersCount = new int [base][base_power_k][base_power_num_others];
		pastOthersCount = new int[base_power_k][base_power_num_others];
		
		// Create constants for tracking prevValues
		maxShiftedValue = new int[base];
		for (int v = 0; v < base; v++) {
			maxShiftedValue[v] = v * MathsUtils.power(base, k-1);
		}
	}

	@Override
	public void initialise(){
		super.initialise();
		
		MatrixUtils.fill(sourceDestPastOthersCount, 0);
		MatrixUtils.fill(sourcePastOthersCount, 0);
		MatrixUtils.fill(destPastOthersCount, 0);
		MatrixUtils.fill(pastOthersCount, 0);
	}
	
	/**
 	 * Add observations for source-destination-conditionals time-series
 	 *  to our estimates of the pdfs.
 	 *  
	 * @param source source time series
	 * @param dest destination time series. Must be of same length as
	 *   source.
	 * @param conditionals conditionals multivariate time series,
	 *  indexed first by time, then by variable number.
	 *  Must be of same length in time as source, and must be
	 *  {@link #numOtherInfoContributors} conditionals here.
	 */
	public void addObservations(int[] source, int[] dest, int[][] conditionals)
		throws Exception {

		int rows = dest.length;
		
		if ((source.length != rows) || (conditionals.length != rows)) {
			throw new Exception("Number of observations must match for dest, source and conditionals");
		}
		
		if (rows - startObservationTime <= 0) {
			return;
		}
		
		// increment the count of observations:
		observations += (rows - startObservationTime); 
		
		if (numOtherInfoContributors != conditionals[0].length) {
			throw new Exception(String.format("conditionals does not have the expected number of variables (%d)", numOtherInfoContributors));
		}

		// Initialise and store the current previous value
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += dest[p];
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal, othersVal;
		for (int r = startObservationTime; r < rows; r++) {
			// Add to the count for this particular transition:
			destVal = dest[r];
			sourceVal = source[r-1];
			othersVal = 0;
			for (int o = 0; o < numOtherInfoContributors; o++) {
				// Include this other contributor
				othersVal *= base;
				othersVal += conditionals[r-1][o];
			}
			sourceDestPastOthersCount[sourceVal][destVal][pastVal][othersVal]++;
			sourcePastOthersCount[sourceVal][pastVal][othersVal]++;
			destPastOthersCount[destVal][pastVal][othersVal]++;
			pastOthersCount[pastVal][othersVal]++;
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[dest[r-k]];
				pastVal *= base;
				pastVal += dest[r];
			}
		}
	}

	/**
 	 * <p>Add observations for a single source-destination-conditionals set
 	 *  to our estimates of the pdfs.</p>
 	 *  
 	 * <p>This method takes in a univariate conditionals time series - it is assumed
 	 * that either numOtherInfoContributors == 1 or the user has combined the 
 	 * multivariate tuples into a single value for each observation, e.g. by calling
 	 * {@link MatrixUtils#computeCombinedValues(int[][], int)}. This cannot be
 	 * checked here however, so use at your own risk!
 	 * </p>
 	 * 
	 * @param source source time series
	 * @param dest destination time series. Must be of same length as
	 *   source.
	 * @param conditionals conditionals univariate time series,
	 *  indexed first by time, then by variable number.
	 *  Must be of same length in time as source, and we must have
	 *  either {@link #numOtherInfoContributors}=1, or the user has
	 *  combined the values of the multivariate conditionals time series
	 *  into single (unique) values at each time step (this is not checked however).
	 */
	public void addObservations(int[] source, int[] dest, int[] conditionals)
		throws Exception {

		int rows = dest.length;
		
		if ((source.length != rows) || (conditionals.length != rows)) {
			throw new Exception("Number of observations must match for dest, source and conditionals");
		}
		
		if (rows - startObservationTime <= 0) {
			return;
		}
		
		// increment the count of observations:
		observations += (rows - startObservationTime); 
		
		// Initialise and store the current previous value
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += dest[p];
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal, othersVal;
		for (int r = startObservationTime; r < rows; r++) {
			// Add to the count for this particular transition:
			destVal = dest[r];
			sourceVal = source[r-1];
			othersVal = conditionals[r-1];
			sourceDestPastOthersCount[sourceVal][destVal][pastVal][othersVal]++;
			sourcePastOthersCount[sourceVal][pastVal][othersVal]++;
			destPastOthersCount[destVal][pastVal][othersVal]++;
			pastOthersCount[pastVal][othersVal]++;
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[dest[r-k]];
				pastVal *= base;
				pastVal += dest[r];
			}
		}
	}

	/**
 	 * Add observations in to our estimates of the pdfs
 	 * from a multivariate time-series of homogeneous variables.
 	 * This call suitable only for homogeneous agents, as all
 	 *  agents will contribute to single pdfs, and all are assumed
 	 *  to have other info contributors at same offsets.
 	 * 
	 * @param states multivariate time series, indexed first by time
	 *  then by variable number.
	 * @param j number of columns to compute transfer entropy across
	 * 	(i.e. for each destination variable i, we have a source
	 *  at i-j, dest i: transfer is j cells to the right) 
	 * @param otherSourcesToDestOffsets offsets of the other information contributors
	 *        from each destination.
	 *        (i.e. offsets from each other information source to the destination -
	 *        the offset is signed the same way as j!)
	 *        othersOffsets is permitted to include j, it will be ignored.
	 */
	public void addObservations(int states[][], int j, int otherSourcesToDestOffsets[]) {
		addObservations(states, j, otherSourcesToDestOffsets, false);
	}
	/**
	 * Private method to implement {@link #addObservations(int[][], int, int[])}
	 * 
	 * @param states multivariate time series, indexed first by time
	 *  then by variable number.
	 * @param j number of columns to compute transfer entropy across
	 * 	(i.e. for each destination variable i, we have a source
	 *  at i-j, dest i: transfer is j cells to the right) 
	 * @param otherSourcesToDestOffsets offsets of the other information contributors
	 *        from each destination.
	 *        (i.e. offsets from each other information source to the destination -
	 *        the offset is signed the same way as j!)
	 *        othersOffsets is permitted to include j, it will be ignored.
	 * @param cleanedOthers whether it has been checked if j is included in otherSourcesToDestOffsets
	 *  or not
	 */
	private void addObservations(int states[][], int j, int otherSourcesToDestOffsets[], boolean cleanedOthers) {
		
		int[] cleanedOthersOffsets;
		if (cleanedOthers) {
			cleanedOthersOffsets = otherSourcesToDestOffsets;
		} else {
			cleanedOthersOffsets = cleanOffsetOthers(otherSourcesToDestOffsets, j, k > 0);
			// This call made redundant by cleanOffsetOthers:
			// confirmEnoughOffsetOthers(othersOffsets, j);
		}
		
		int rows = states.length;
		int columns = states[0].length;
		// increment the count of observations:
		observations += (rows - startObservationTime)*columns; 
		
		// Initialise and store the current previous value for each column
		int[] pastVal = new int[columns]; 
		for (int c = 0; c < columns; c++) {
			pastVal[c] = 0;
			for (int p = 0; p < k; p++) {
				pastVal[c] *= base;
				pastVal[c] += states[p][c];
			}
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal, othersVal;
		for (int r = startObservationTime; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				// Add to the count for this particular transition:
				// (cell's assigned as above)
				destVal = states[r][c];
				sourceVal = states[r-1][(c-j+columns) % columns];
				othersVal = 0;
				for (int o = 0; o < cleanedOthersOffsets.length; o++) {
					// Include this other contributor
					othersVal *= base;
					othersVal += states[r-1][(c-cleanedOthersOffsets[o]+columns) % columns];
				}
				sourceDestPastOthersCount[sourceVal][destVal][pastVal[c]][othersVal]++;
				sourcePastOthersCount[sourceVal][pastVal[c]][othersVal]++;
				destPastOthersCount[destVal][pastVal[c]][othersVal]++;
				pastOthersCount[pastVal[c]][othersVal]++;
				// Update the previous value:
				if (k > 0) {
					pastVal[c] -= maxShiftedValue[states[r-k][c]];
					pastVal[c] *= base;
					pastVal[c] += states[r][c];
				}
			}
		}		
	}
	
	/**
 	 * Add observations for a single source-destination-conditionals set of the
 	 *  multivariate time series.
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to addObservations(int states[][])
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states multivariate time series, indexed first by time
	 *  then by variable number.
	 * @param sourceCol column index of the source
	 * @param destCol column index of the destination
	 * @param othersAbsolute column indices of the conditional variables.
	 *  othersAbsolute is permitted to include sourceCol or destCol (if k>0),
	 *  they will be ignored.
	 */
	public void addObservations(int states[][], int sourceCol, int destCol, int[] othersAbsolute) {
		addObservations(states, sourceCol, destCol, othersAbsolute, false);
	}
	/**
	 * Private method to implement {@link #addObservations(int[][], int, int, int[])}
	 * 
	 * @param states multivariate time series, indexed first by time
	 *  then by variable number.
	 * @param sourceCol column index of the source
	 * @param destCol column index of the destination
	 * @param othersAbsolute column indices of the conditional variables.
	 *  othersAbsolute is permitted to include sourceCol or destCol (if k>0),
	 *  they will be ignored.
	 * @param cleanedOthers whether it has been checked if othersAbsolute
	 *  contains sourceCol or destCol
	 */
	private void addObservations(int states[][], int sourceCol, int destCol, int[] othersAbsolute, boolean cleanedOthers) {

		int[] cleanedOthersAbsolute;
		if (cleanedOthers) {
			cleanedOthersAbsolute = othersAbsolute;
		} else {
			cleanedOthersAbsolute = cleanAbsoluteOthers(othersAbsolute, sourceCol,
					destCol, k > 0);
			// This call made redundant by cleanAbsoluteOthers:
			// confirmEnoughAbsoluteOthers(othersAbsolute, destCol, sourceCol);
		}
		
		int rows = states.length;
		// increment the count of observations:
		observations += (rows - startObservationTime); 
		
		// Initialise and store the current previous value for each column
		int pastVal = 0;
		for (int p = 0; p < k; p++) {
			pastVal *= base;
			pastVal += states[p][destCol];
		}
		
		// 1. Count the tuples observed
		int destVal, sourceVal, othersVal;
		for (int r = startObservationTime; r < rows; r++) {
			// Add to the count for this particular transition:
			// (cell's assigned as above)
			destVal = states[r][destCol];
			sourceVal = states[r-1][sourceCol];
			othersVal = 0;
			for (int o = 0; o < cleanedOthersAbsolute.length; o++) {
				// Include this other contributor
				othersVal *= base;
				othersVal += states[r-1][cleanedOthersAbsolute[o]];
			}
			sourceDestPastOthersCount[sourceVal][destVal][pastVal][othersVal]++;
			sourcePastOthersCount[sourceVal][pastVal][othersVal]++;
			destPastOthersCount[destVal][pastVal][othersVal]++;
			pastOthersCount[pastVal][othersVal]++;
			// Update the previous value:
			if (k > 0) {
				pastVal -= maxShiftedValue[states[r-k][destCol]];
				pastVal *= base;
				pastVal += states[r][destCol];
			}
		}
	}

	@Override
	public double computeAverageLocalOfObservations() {
		double te = 0.0;
		double teCont = 0.0;

		max = 0;
		min = 0;
		double meanSqLocals = 0;
		for (int othersVal = 0; othersVal < this.base_power_num_others; othersVal++) {
			for (int pastVal = 0; pastVal < base_power_k; pastVal++) {
				for (int destVal = 0; destVal < base; destVal++) {
					for (int sourceVal = 0; sourceVal < base; sourceVal++) {
						// Compute TE contribution:
						if (sourceDestPastOthersCount[sourceVal][destVal][pastVal][othersVal] != 0) {
							/* Double check: should never happen
							if ((sourcePastCount[sourceVal][pastVal][othersVal] == 0) ||
								(destPastCount[destVal][pastVal][othersVal] == 0) ||
								(pastCount[pastVal][othersVal] == 0)) {
								throw new RuntimeException("one subcount was zero!!");
							}
							*/
							// compute p(source,dest,past)
							double p_source_dest_past_others = (double)
								sourceDestPastOthersCount[sourceVal][destVal][pastVal][othersVal] / (double) observations;
							
							double logTerm = ((double) sourceDestPastOthersCount[sourceVal][destVal][pastVal][othersVal] / (double) sourcePastOthersCount[sourceVal][pastVal][othersVal]) /
							 	((double) destPastOthersCount[destVal][pastVal][othersVal] / (double) pastOthersCount[pastVal][othersVal]);
							double localValue = Math.log(logTerm) / log_2;
							teCont = p_source_dest_past_others * localValue;
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
		}
		
		average = te;
		std = Math.sqrt(meanSqLocals - average * average);
		return te;
	}
	
	/**
	 * Generate a bootstrapped distribution of what the 
	 * conditional TE would look like,
	 * under a null hypothesis that the source values of our
	 * samples had no relation to the destination value
	 * (in the context of the destination past and conditionals).
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for MI,
	 * conditional MI and TE.
	 * </p>
	 * 
	 * <p>Note that if several disjoint time-series have been added 
	 * as observations using {@link #addObservations(int[], int[])} etc.,
	 * then these separate "trials" will be mixed up in the generation
	 * of surrogates here.</p>
	 * 
	 * @param numPermutationsToCheck number of surrogate samples to bootstrap
	 *  to generate the distribution.
	 * @return the distribution of conditional TE scores under this null hypothesis.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) {
		double actualTE = computeAverageLocalOfObservations();
		
		// Reconstruct the source values (not necessarily in order)
		int[] sourceValues = new int[observations];
		int t_s = 0;
		for (int sourceVal = 0; sourceVal < base; sourceVal++) {
			// Count up the number of times this source value was observed:
			int numberOfSamples = 0;
			for (int pastVal = 0; pastVal < base_power_k; pastVal++) {
				for (int othersVal = 0; othersVal < this.base_power_num_others; othersVal++) {
					numberOfSamples += sourcePastOthersCount[sourceVal][pastVal][othersVal];
				}
			}
			// Now add all of these as unordered observations:
			MatrixUtils.fill(sourceValues, sourceVal, t_s, numberOfSamples);
			t_s += numberOfSamples;
		}
		
		// And construct unordered (dest,past,others) tuples.
		// It doesn't matter that we've appeared to destroy the ordering here because
		//  the joint distribution destPastOthersCount is actually preserved in
		//  our construction of pastVal and destValues and othersValues together here.
		int[] destValues = new int[observations];
		int[] pastValues = new int[observations];
		int[] othersValues = new int[observations];
		int t_d = 0;
		int t_p = 0;
		int t_o = 0;
		for (int pastVal = 0; pastVal < base_power_k; pastVal++) {
			for (int othersVal = 0; othersVal < this.base_power_num_others; othersVal++) {
				// Add in pastOthersCount[pastVal][othersVal] dummy past values
				MatrixUtils.fill(pastValues, pastVal, t_p,
						pastOthersCount[pastVal][othersVal]);
				t_p += pastOthersCount[pastVal][othersVal];
				// Add in pastOthersCount[pastVal][othersVal] dummy others values
				MatrixUtils.fill(othersValues, othersVal, t_o,
						pastOthersCount[pastVal][othersVal]);
				t_o += pastOthersCount[pastVal][othersVal];
				for (int destVal = 0; destVal < base; destVal++) {
					MatrixUtils.fill(destValues, destVal, t_d,
						destPastOthersCount[destVal][pastVal][othersVal]);
					t_d += destPastOthersCount[destVal][pastVal][othersVal];
				}
			}
		}
		
		// Construct new source orderings based on the source probabilities only
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(observations, numPermutationsToCheck);

		ConditionalTransferEntropyCalculatorDiscrete cte = newInstance(base, k, numOtherInfoContributors);
		cte.initialise();
		cte.observations = observations;
		cte.pastOthersCount = pastOthersCount;
		cte.destPastOthersCount = destPastOthersCount;
		int countWhereTeIsMoreSignificantThanOriginal = 0;
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		for (int p = 0; p < numPermutationsToCheck; p++) {
			// Generate a new re-ordered data set for the source
			int[] newSourceData = MatrixUtils.extractSelectedTimePoints(sourceValues, newOrderings[p]);
			// compute the joint probability distributions
			MatrixUtils.fill(cte.sourceDestPastOthersCount, 0);
			MatrixUtils.fill(cte.sourcePastOthersCount, 0);
			for (int t = 0; t < observations; t++) {
				cte.sourcePastOthersCount[newSourceData[t]][pastValues[t]][othersValues[t]]++;
				cte.sourceDestPastOthersCount[newSourceData[t]][destValues[t]]
				                            [pastValues[t]][othersValues[t]]++;
			}
			// And get a TE value for this realisation:
			double newTe = cte.computeAverageLocalOfObservations();
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
	
	/**
	 * Computes local conditional transfer entropy for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method to be used for homogeneous agents only
	 *  
	 * @param states multivariate time series, indexed first by time
	 *  then by variable number.
	 * @param j number of columns to compute transfer entropy across
	 * 	(i.e. for each destination variable i, we have a source
	 *  at i-j, dest i: transfer is j cells to the right) 
	 * @param otherSourcesToDestOffsets offsets of the other information contributors
	 *        from each destination.
	 *        (i.e. offsets from each other information source to the destination -
	 *        the offset is signed the same way as j!)
	 *        othersOffsets is permitted to include j, it will be ignored.
	 * @return 2D time-series of local conditional TE values (indexed 
	 *  as per states)
	 */
	public double[][] computeLocalFromPreviousObservations
		(int states[][], int j, int otherSourcesToDestOffsets[]){
		
		return computeLocalFromPreviousObservations(states, j, otherSourcesToDestOffsets, false);
	}
	/**
	 * Private method to implement {@link #computeLocalFromPreviousObservations(int[][], int, int[])}
	 * 
	 * @param states multivariate time series, indexed first by time
	 *  then by variable number.
	 * @param j number of columns to compute transfer entropy across
	 * 	(i.e. for each destination variable i, we have a source
	 *  at i-j, dest i: transfer is j cells to the right) 
	 * @param otherSourcesToDestOffsets offsets of the other information contributors
	 *        from each destination.
	 *        (i.e. offsets from each other information source to the destination -
	 *        the offset is signed the same way as j!)
	 *        othersOffsets is permitted to include j, it will be ignored.
	 * @param cleanedOthers whether it has been checked if j is included in otherSourcesToDestOffsets
	 *  or not
	 * @return 2D time-series of local conditional TE values (indexed 
	 *  as per states)
	 */
	private double[][] computeLocalFromPreviousObservations
		(int states[][], int j, int othersOffsets[], boolean cleanedOthers){
		
		int[] cleanedOthersOffsets;
		if (cleanedOthers) {
			cleanedOthersOffsets = othersOffsets;
		} else {
			cleanedOthersOffsets = cleanOffsetOthers(othersOffsets, j, k > 0);
			// This call made redundant by cleanOffsetOthers:
			// confirmEnoughOffsetOthers(othersOffsets, j);
		}

		int rows = states.length;
		int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[][] localTE = new double[rows][columns];
		average = 0;
		max = 0;
		min = 0;

		// Initialise and store the current previous value for each column
		int[] pastVal = new int[columns]; 
		for (int c = 0; c < columns; c++) {
			pastVal[c] = 0;
			for (int p = 0; p < k; p++) {
				pastVal[c] *= base;
				pastVal[c] += states[p][c];
			}
		}
		int destVal, sourceVal, othersVal;
		double logTerm;
		for (int r = startObservationTime; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				destVal = states[r][c];
				sourceVal = states[r-1][(c-j+columns) % columns];
				othersVal = 0;
				for (int o = 0; o < cleanedOthersOffsets.length; o++) {
					// Include this other contributor
					othersVal *= base;
					othersVal += states[r-1][(c-cleanedOthersOffsets[o]+columns) % columns];
				}

				// Now compute the local value
				logTerm = ((double) sourceDestPastOthersCount[sourceVal][destVal][pastVal[c]][othersVal] / (double) sourcePastOthersCount[sourceVal][pastVal[c]][othersVal]) /
			 		((double) destPastOthersCount[destVal][pastVal[c]][othersVal] / (double) pastOthersCount[pastVal[c]][othersVal]);
				localTE[r][c] = Math.log(logTerm) / log_2;
				average += localTE[r][c];
				if (localTE[r][c] > max) {
					max = localTE[r][c];
				} else if (localTE[r][c] < min) {
					min = localTE[r][c];
				}
				// Update the previous value:
				if (k > 0) {
					pastVal[c] -= maxShiftedValue[states[r-k][c]];
					pastVal[c] *= base;
					pastVal[c] += states[r][c];
				}
			}
		}		

		average = average/(double) (columns * (rows - startObservationTime));
		
		return localTE;
	}
	
	/**
	 * Computes local active information storage for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method is suitable for heterogeneous agents
	 *  
	 * @param states multivariate time series, indexed first by time
	 *  then by variable number.
	 * @param sourceCol column index of the source
	 * @param destCol column index of the destination
	 * @param othersAbsolute column indices of the conditional variables.
	 *  othersAbsolute is permitted to include sourceCol or destCol (if k>0),
	 *  they will be ignored.
	 * @return time-series of local conditional TE values
	 */
	public double[] computeLocalFromPreviousObservations
		(int states[][], int sourceCol, int destCol, int[] othersAbsolute){
		
		return computeLocalFromPreviousObservations(states, sourceCol, destCol, othersAbsolute, false);
	}
	/**
	 * Private method to implement {@link #computeLocalFromPreviousObservations(int[][], int, int, int[])}
	 * 
	 * @param states multivariate time series, indexed first by time
	 *  then by variable number.
	 * @param sourceCol column index of the source
	 * @param destCol column index of the destination
	 * @param othersAbsolute column indices of the conditional variables.
	 *  othersAbsolute is permitted to include sourceCol or destCol (if k>0),
	 *  they will be ignored.
	 * @param cleanedOthers whether it has been checked if othersAbsolute
	 *  contains sourceCol or destCol
	 * @return time-series of local conditional TE values
	 */
	private double[] computeLocalFromPreviousObservations
		(int states[][], int sourceCol, int destCol, int[] othersAbsolute, boolean cleanedOthers){

		int[] cleanedOthersAbsolute;
		if (cleanedOthers) {
			cleanedOthersAbsolute = othersAbsolute;
		} else {
			cleanedOthersAbsolute = cleanAbsoluteOthers(othersAbsolute,
					sourceCol, destCol, k > 0);
			// This call made redundant by cleanOffsetOthers:
			// confirmEnoughAbsoluteOthers(othersAbsolute, destCol, sourceCol);
		}
		
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
			pastVal += states[p][destCol];
		}
		int destVal, sourceVal, othersVal;
		double logTerm;
		for (int r = startObservationTime; r < rows; r++) {
			destVal = states[r][destCol];
			sourceVal = states[r-1][sourceCol];
			othersVal = 0;
			for (int o = 0; o < cleanedOthersAbsolute.length; o++) {
				// Include this other contributor
				othersVal *= base;
				othersVal += states[r-1][cleanedOthersAbsolute[o]];
			}
			// Now compute the local value
			logTerm = ((double) sourceDestPastOthersCount[sourceVal][destVal][pastVal][othersVal] / (double) sourcePastOthersCount[sourceVal][pastVal][othersVal]) /
		 		((double) destPastOthersCount[destVal][pastVal][othersVal] / (double) pastOthersCount[pastVal][othersVal]);
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
		}

		average = average/(double) (rows - startObservationTime);
		
		return localTE;
	}

	/**
	 * Standalone routine to 
	 * compute local conditional transfer entropy across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * This method to be called for homogeneous agents only
	 * 
	 * @param states multivariate time series, indexed first by time
	 *  then by variable number.
	 * @param j number of columns to compute transfer entropy across
	 * 	(i.e. for each destination variable i, we have a source
	 *  at i-j, dest i: transfer is j cells to the right) 
	 * @param otherSourcesToDestOffsets offsets of the other information contributors
	 *        from each destination.
	 *        (i.e. offsets from each other information source to the destination -
	 *        the offset is signed the same way as j!)
	 *        othersOffsets is permitted to include j, it will be ignored.
	 * @return 2D time-series of local conditional TE values (indexed 
	 *  as per states)
	 */
	public double[][] computeLocal(int states[][], int j, int[] otherSourcesToDestOffsets) {
		
		initialise();
		int[] cleanedOthersOffsets = cleanOffsetOthers(otherSourcesToDestOffsets, j, k > 0);
		addObservations(states, j, cleanedOthersOffsets, true);
		return computeLocalFromPreviousObservations(states, j, cleanedOthersOffsets, true);
	}
	
	/**
	 * Standalone routine to 
	 * compute average conditional transfer entropy across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return the average
	 * This method to be called for homogeneous agents only
	 * 
	 * @param states multivariate time series, indexed first by time
	 *  then by variable number.
	 * @param j number of columns to compute transfer entropy across
	 * 	(i.e. for each destination variable i, we have a source
	 *  at i-j, dest i: transfer is j cells to the right) 
	 * @param otherSourcesToDestOffsets offsets of the other information contributors
	 *        from each destination.
	 *        (i.e. offsets from each other information source to the destination -
	 *        the offset is signed the same way as j!)
	 *        othersOffsets is permitted to include j, it will be ignored.
	 * @return average conditional TE from these observations
	 */
	public double computeAverageLocal(int states[][], int j, int[] otherSourcesToDestOffsets) {
		
		initialise();
		addObservations(states, j, otherSourcesToDestOffsets);
		return computeAverageLocalOfObservations();
	}

	/**
	 * Standalone routine to 
	 * compute local conditional transfer entropy for a specific set of
	 * source-destination-conditionals in a multivariate time-series.
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros.
	 * This method suitable for heterogeneous agents.
	 * 
	 * @param states multivariate time series, indexed first by time
	 *  then by variable number.
	 * @param sourceCol column index of the source
	 * @param destCol column index of the destination
	 * @param othersAbsolute column indices of the conditional variables.
	 *  othersAbsolute is permitted to include sourceCol or destCol (if k>0),
	 *  they will be ignored.
	 * @return time-series of local conditional TE values for these 
	 *  observations
	 */
	public double[] computeLocal(int states[][], int sourceCol, int destCol, int[] othersAbsolute) {
		
		initialise();
		int[] cleanedOthers = cleanAbsoluteOthers(othersAbsolute, sourceCol,
				destCol, k > 0);
		addObservations(states, sourceCol, destCol, cleanedOthers, true);
		return computeLocalFromPreviousObservations(states, sourceCol, destCol, cleanedOthers, true);
	}

	/**
	 * Standalone routine to 
	 * compute average conditional transfer entropy for a specific set of
	 * source-destination-conditionals in a multivariate time-series.
	 * This method suitable for heterogeneous agents
	 * 
	 * @param states multivariate time series, indexed first by time
	 *  then by variable number.
	 * @param sourceCol column index of the source
	 * @param destCol column index of the destination
	 * @param othersAbsolute column indices of the conditional variables.
	 *  othersAbsolute is permitted to include sourceCol or destCol (if k>0),
	 *  they will be ignored.
	 * @return average conditional transfer entropy for these observations
	 */
	public double computeAverageLocal(int states[][], int sourceCol, int destCol, int[] othersAbsolute) {
		
		initialise();
		addObservations(states, sourceCol, destCol, othersAbsolute);
		return computeAverageLocalOfObservations();
	}
	
	/**
	 * Counts the unique information contributors to this node which
	 * are not equal to the source at offset j from the destination,
	 * or the node itself (offset 0,
	 * node itself not included only when removeDest is set to true).
	 * 
	 * <p>This is primarily intended for use inside the method, but
	 * made public as a utility.</p>
	 * 
	 * @param otherSourcesToDestOffsets array of offsets of the destination from each source
	 * @param j offset of the destination from the source
	 * @param removeDest remove the destination itself from the count
	 *     of offset others.
	 * @return the number of unique such sources
	 */
	public static int countOfOffsetOthers(int[] otherSourcesToDestOffsets, int j,
			boolean removeDest) {
		int countOfOthers = 0;
		for (int index = 0; index < otherSourcesToDestOffsets.length; index++) {
			if ((otherSourcesToDestOffsets[index] != j) && 
				((otherSourcesToDestOffsets[index] != 0) || !removeDest)) {
				countOfOthers++;
			}
		}
		return countOfOthers;
	}
	
	/**
	 * Counts the information contributors to the dest which
	 * are not equal to src or the node itself (offset 0,
	 * node itself not included only when removeDest is set to true)
	 * 
	 * <p>This is primarily intended for use inside the method, but
	 * made public as a utility.</p>
	 * 
	 * @param others array of source indices
	 * @param src current source index we are considering
	 * @param dest current dest index we are considering
	 * @param removeDest remove the destination itself from the count
	 *     of absolute others.
	 * @return the number of unique such sources
	 */
	public static int countOfAbsoluteOthers(int[] others, int src, int dest,
			boolean removeDest) {
		int countOfOthers = 0;
		for (int index = 0; index < others.length; index++) {
			if ((others[index] != src) &&
				((others[index] != dest) || !removeDest)) {
				countOfOthers++;
			}
		}
		return countOfOthers;
	}

	/**
	 * Check that the supplied array of offsets as other info 
	 * contributors is long enough compared to our expectation
	 * 
	 * @param othersOffsets array of offsets from each source to the destination
	 * @param j offset from the source to the destination
	 * @param removeDest remove the destination itself from the count
	 *     of absolute others.
	 * @return whether the number of such sources matches what we expect
	 * @throws Exception if the number of such sources does not match.
	 */
	public boolean confirmEnoughOffsetOthers(int[] othersOffsets, int j,
			boolean removeDest) {
		if (countOfOffsetOthers(othersOffsets, j, removeDest) !=
				numOtherInfoContributors) {
			throw new RuntimeException("Incorrect number of others in offsets");
		}
		return true;
	}

	/**
	 * Check that the supplied array of absolutes as other info 
	 * contributors is long enough compared to our expectation
	 * 
	 * @param othersAbsolute array of source indices
	 * @param src current source index we are considering
	 * @param dest current dest index we are considering
	 * @param removeDest remove the destination itself from the count
	 *     of absolute others.
	 * @return whether the number of such sources matches what we expect
	 * @throws Exception if the number of such sources does not match.
	 */
	public boolean confirmEnoughAbsoluteOthers(int[] othersAbsolute, int src,
			int dest, boolean removeDest) {
		if (countOfAbsoluteOthers(othersAbsolute, src, dest, removeDest) !=
				numOtherInfoContributors) {
			throw new RuntimeException("Incorrect number of others in absolutes");
		}
		return true;
	}
	
	/**
	 * Returns the information contributors to this node which
	 * are not equal to the offset j or the node itself (offset 0,
	 * removed only if removeDest is set to true).
	 * Checks that there are enough other information contributors.
	 * 
	 * @param othersOffsets array of offsets from each source to the destination
	 * @param j offset from the source to the destination
	 * @param removeDest Remove the destination itself from the cleaned
	 *           other sources (if it is there). Should not be done
	 *           if k == 0 (because then the destination is not included
	 *           in the past history)
	 * @return othersOffsets with entries for the offending 
	 *  contributors removed (i.e. array may be shortened)
	 */
	public int[] cleanOffsetOthers(int[] othersOffsets, int j, boolean removeDest) {
		int[] cleaned = new int[numOtherInfoContributors];
		int countOfOthers = 0;
		for (int index = 0; index < othersOffsets.length; index++) {
			if ((othersOffsets[index] != j) &&
				((othersOffsets[index] != 0) || !removeDest)) {
				// Add this candidate source to the cleaned sources
				if (countOfOthers == numOtherInfoContributors) {
					// We've already taken all the other info 
					//  contributors we expected
					countOfOthers++;
					break;
				}
				cleaned[countOfOthers] = othersOffsets[index];
				countOfOthers++;
			}
		}
		if (countOfOthers < numOtherInfoContributors) {
			throw new RuntimeException("Too few others in offsets");
		} else if (countOfOthers > numOtherInfoContributors) {
			throw new RuntimeException("Too many others in offsets");
		}
		return cleaned;
	}
	
	/**
	 * Returns the information contributors to the dest which
	 * are not equal to src or the node itself (offset 0,
	 * removed only if removeDest is true).
	 * Checks that there are enough other information contributors.
	 * 
	 * @param others array of source indices
	 * @param src current source index we are considering
	 * @param dest current dest index we are considering
	 * @param removeDest Remove the destination itself from the cleaned
	 *           other sources (if it is there). Should not be done
	 *           if k == 0 (because then the destination is not included
	 *           in the past history)
	 * @return others with entries for the offending 
	 *  contributors removed (i.e. array may be shortened)
	 */
	public int[] cleanAbsoluteOthers(int[] others, int src, int dest,
			boolean removeDest) {
		int[] cleaned = new int[numOtherInfoContributors];
		int countOfOthers = 0;
		for (int index = 0; index < others.length; index++) {
			if ((others[index] != src) &&
				((others[index] != dest) || !removeDest)) {
				// Add this candidate source to the cleaned sources
				if (countOfOthers == numOtherInfoContributors) {
					// We've already taken all the other info 
					//  contributors we expected
					countOfOthers++;
					break;
				}
				cleaned[countOfOthers] = others[index];
				countOfOthers++;
			}
		}
		if (countOfOthers < numOtherInfoContributors) {
			throw new RuntimeException("Too few others in absolutes");
		} else if (countOfOthers > numOtherInfoContributors) {
			throw new RuntimeException("Too many others in absolutes");
		}
		return cleaned;
	}
}
