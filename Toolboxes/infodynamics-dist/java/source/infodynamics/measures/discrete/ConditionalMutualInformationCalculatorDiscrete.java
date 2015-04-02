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
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

/**
 * <p>Conditional Mutual information (MI) calculator for univariate discrete (int[]) data.</p>
 * 
 * <p>Usage of the class is intended to follow this paradigm:</p>
 * <ol>
 * 		<li>Construct the calculator:
 * 			{@link #ConditionalMutualInformationCalculatorDiscrete(int, int, int)};</li>
 *		<li>Initialise the calculator using {@link #initialise()};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using one or more calls to
 * 			sets of {@link #addObservations(int[], int[], int[])} methods, then</li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average MI: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>local MI values, such as
 * 				{@link #computeLocalFromPreviousObservations(int[], int[], int[])};</li>
 * 				<li>comparison to null distribution, such as
 * 				{@link #computeSignificance()};</li>
 * 				<li>and variants of these.</li>
 * 			</ul>
 * 		</li>
 * 		<li>As an alternative to steps 3 and 4, the user may undertake
 * 			standalone computation from a single set of observations, via
 *  		e.g.: {@link #computeLocal(int[], int[], int[])}.</li>
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
public class ConditionalMutualInformationCalculatorDiscrete
	extends InfoMeasureCalculatorDiscrete implements AnalyticNullDistributionComputer {

	/**
	 * Store the number of symbols for each variable
	 */
	protected int base1;
	protected int base2;
	protected int condBase;
	
	protected int[][][] firstSecondCondCount = null;	// count for (x,y,Cond) tuples
	protected int[][] firstCondCount = null;			// count for (x,Cond) tuples
	protected int[][] secondCondCount = null; // Count for (y,Cond) tuples
	protected int[] condCount = null; // Count for Cond

	protected boolean condMiComputed = false;
	
	/**
	 * User was formerly forced to create new instances through this factory method.
	 * Retained for backwards compatibility.
	 * 
	 * @param base1 number of symbols for first variable.
	 *        E.g. binary variables are in base-2.
	 * @param base2 number of symbols for second variable.
	 * @param condBase number of symbols for conditional variable.
	 * @deprecated As of JIDT 1.0, call {@link #ConditionalMutualInformationCalculator(int, int, int)}
	 * 				directly.
	 * @return new calculator object
	 */
	public static ConditionalMutualInformationCalculatorDiscrete newInstance(int base1, int base2, int condBase) {
		return new ConditionalMutualInformationCalculatorDiscrete(base1, base2, condBase);
	}

	/**
	 * Construct a new instance
	 * 
	 * @param base1 number of symbols for first variable.
	 *        E.g. binary variables are in base-2.
	 * @param base2 number of symbols for second variable.
	 * @param condBase number of symbols for conditional variable.
	 */
	public ConditionalMutualInformationCalculatorDiscrete(int base1, int base2, int condBase) {

		// Create super object, just with first base
		super(base1);
		
		// Store the bases
		this.base1 = base1;
		this.base2 = base2;
		this.condBase = condBase;
		
		// Create storage for extra counts of observations
		firstSecondCondCount = new int[base1][base2][condBase];
		firstCondCount = new int[base1][condBase];		
		secondCondCount = new int[base2][condBase];
		condCount = new int[condBase];
	}

	@Override
	public void initialise(){
		super.initialise();
		
		condMiComputed = false;
		
		MatrixUtils.fill(firstSecondCondCount, 0);
		MatrixUtils.fill(firstCondCount, 0);
		MatrixUtils.fill(secondCondCount, 0);
		MatrixUtils.fill(condCount,0);
	}
		
	/**
 	 * Add observations for the given var1,var2,cond tuples
 	 *  of the variables
 	 *  to our estimates of the pdfs.
	 *
	 * @param var1 series of values for the first variable
	 * @param var2 series of values for the second variable
	 * @param cond series of values for the conditional variable
	 */
	public void addObservations(int var1[], int var2[], int cond[]) {
		int rows = var1.length;
		// increment the count of observations:
		observations += rows;
		
		// 1. Count the tuples observed
		for (int r = 0; r < rows; r++) {
			// Add to the count for this particular transition:
			firstSecondCondCount[var1[r]][var2[r]][cond[r]]++;
			firstCondCount[var1[r]][cond[r]]++;
			secondCondCount[var2[r]][cond[r]]++;
			condCount[cond[r]]++;
		}
	}

	/**
 	 * Add observations for the given var1,var2,cond tuples
 	 *  of the variables
 	 *  to our estimates of the pdfs.
 	 * This method signature allows applications using byte data (to save memory)
 	 *  to operate correctly.
	 *
	 * @param var1 series of values for the first variable
	 * @param var2 series of values for the second variable
	 * @param cond series of values for the conditional variable
	 */
	public void addObservations(byte var1[], byte var2[], byte cond[]) {
		int rows = var1.length;
		// increment the count of observations:
		observations += rows;
		
		// 1. Count the tuples observed
		for (int r = 0; r < rows; r++) {
			// Add to the count for this particular transition:
			firstSecondCondCount[var1[r]][var2[r]][cond[r]]++;
			firstCondCount[var1[r]][cond[r]]++;
			secondCondCount[var2[r]][cond[r]]++;
			condCount[cond[r]]++;
		}
	}

	/**
 	 * Add observations for the given var1,var2,cond tuples
 	 *  of the variables
 	 *  to our estimates of the pdfs, only when those observations are confirmed
 	 *  as valid.
	 *
	 * @param var1 series of values for the first variable
	 * @param var2 series of values for the second variable
	 * @param cond series of values for the conditional variable
	 * @param valid series of whether each observation is valid
	 * 	to be counted in the PDFs
	 */
	public void addObservations(int var1[], int var2[], int cond[],
			boolean[] valid) {
		
		int rows = var1.length;
		
		// 1. Count the tuples observed
		for (int r = 0; r < rows; r++) {
			if (valid[r]) {
				// Add to the count for this particular transition:
				firstSecondCondCount[var1[r]][var2[r]][cond[r]]++;
				firstCondCount[var1[r]][cond[r]]++;
				secondCondCount[var2[r]][cond[r]]++;
				condCount[cond[r]]++;
				observations++; // increment the count of observations:
			}
		}
	}

	/**
 	 * Add observations for the given var1,var2,cond tuples
 	 *  of the variables
 	 *  to our estimates of the pdfs, only when those observations are confirmed
 	 *  as valid.
 	 * This method signature allows applications using byte data (to save memory)
 	 *  to operate correctly.
	 *
	 * @param var1 series of values for the first variable
	 * @param var2 series of values for the second variable
	 * @param cond series of values for the conditional variable
	 * @param valid whether each observation is valid
	 * 	to be counted in the PDFs
	 */
	public void addObservations(byte var1[], byte var2[], byte cond[],
			boolean[] valid) {
		
		int rows = var1.length;
		
		// 1. Count the tuples observed
		for (int r = 0; r < rows; r++) {
			if (valid[r]) {
				// Add to the count for this particular transition:
				firstSecondCondCount[var1[r]][var2[r]][cond[r]]++;
				firstCondCount[var1[r]][cond[r]]++;
				secondCondCount[var2[r]][cond[r]]++;
				condCount[cond[r]]++;
				observations++; // increment the count of observations:
			}
		}
	}

	/**
 	 * Add observations for the given var1,var2,cond tuples
 	 * 	of the variables
 	 *  to our estimates of the pdfs.
 	 * This method signature allows multiple trials to have their
 	 *  observations added to the PDFs by one method call.
	 *
	 * @param var1 series of values for the first variable;
	 *   first index is time or observation number, second
	 *   index is trial (a second observation number if you will)
	 * @param var2 series of values for the second variable,
	 * 	indexed as per var1
	 * @param cond series of values for the conditional variable,
	 * 	indexed as per var1
	 */
	public void addObservations(int var1[][], int var2[][], int cond[][]) {
		int rows = var1.length;
		int cols = var1[0].length;
		// increment the count of observations:
		observations += rows * cols;
		
		// 1. Count the tuples observed
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < cols; c++) {
				// Add to the count for this particular transition:
				firstSecondCondCount[var1[r][c]][var2[r][c]][cond[r][c]]++;
				firstCondCount[var1[r][c]][cond[r][c]]++;
				secondCondCount[var2[r][c]][cond[r][c]]++;
				condCount[cond[r][c]]++;
			}
		}
	}

	@Override
	public double computeAverageLocalOfObservations() {
		double condMi = 0.0;
		double condMiCont = 0.0;

		max = 0;
		min = 0;
		double meanSqLocals = 0;
		for (int condVal = 0; condVal < condBase; condVal++) {
			// compute p(cond)
			// double p_cond = (double) condCount[condVal] / (double) observations;
			for (int var2Val = 0; var2Val < base2; var2Val++) {
				// compute p(var2,cond)
				// double p_var2_cond = (double) seondCondCount[var2Val][condVal] / (double) observations;
				for (int var1Val = 0; var1Val < base1; var1Val++) {
					// compute p(var1,var2,cond)
					double p_var1_var2_cond = (double) firstSecondCondCount[var1Val][var2Val][condVal] / (double) observations;
					// compute p(var1,cond)
					// double p_var1_cond = (double) firstCondCount[var1Val][condVal] / (double) observations;
					// Compute TE contribution:
					if (firstSecondCondCount[var1Val][var2Val][condVal] != 0) {
						/* Double check: should never happen
						if ((sourcePastCount[sourceVal][pastVal] == 0) ||
							(destPastCount[destVal][pastVal] == 0) ||
							(pastCount[pastVal] == 0)) {
							throw new RuntimeException("one subcount was zero!!");
						}
						*/
						
						double logTerm = ((double) firstSecondCondCount[var1Val][var2Val][condVal] / (double) firstCondCount[var1Val][condVal]) /
						 	((double) secondCondCount[var2Val][condVal] / (double) condCount[condVal]);
						double localValue = Math.log(logTerm) / log_2;
						condMiCont = p_var1_var2_cond * localValue;
						if (localValue > max) {
							max = localValue;
						} else if (localValue < min) {
							min = localValue;
						}
						// Add this contribution to the mean 
						//  of the squared local values
						meanSqLocals += condMiCont * localValue;
					} else {
						condMiCont = 0.0;
					}
					condMi += condMiCont;
				}
			}
		}
		
		average = condMi;
		std = Math.sqrt(meanSqLocals - average * average);
		condMiComputed = true;
		return condMi;
	}
	
	/**
	 * Dump a debug print of the PDFs of our observations
	 */
	public void debugPrintObservations() {
		System.out.println("Var1\tVar2\tCond\tc(1,2,c)\tc(1,c)\tc(2,c)\tc(c)");
		for (int condVal = 0; condVal < condBase; condVal++) {
			// compute p(cond)
			// double p_cond = (double) condCount[condVal] / (double) observations;
			for (int var2Val = 0; var2Val < base2; var2Val++) {
				// compute p(var2,cond)
				// double p_var2_cond = (double) seondCondCount[var2Val][condVal] / (double) observations;
				for (int var1Val = 0; var1Val < base1; var1Val++) {
					// compute p(var1,var2,cond)
					// double p_var1_var2_cond = (double) firstSecondCondCount[var1Val][var2Val][condVal] / (double) observations;
					// compute p(var1,cond)
					// double p_var1_cond = (double) firstCondCount[var1Val][condVal] / (double) observations;
					// Compute TE contribution:
					System.out.println(var1Val + "\t" + var2Val + "\t" + condVal + "\t" +
							firstSecondCondCount[var1Val][var2Val][condVal] + "\t\t" +
							firstCondCount[var1Val][condVal] + "\t" + 
							secondCondCount[var2Val][condVal] + "\t" +
							condCount[condVal]);
				}
			}
		}
	}

	/**
	 * Computes local conditional MI for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 *  
	 * @param var1 series of values for the first variable
	 * @param var2 series of values for the second variable
	 * @param cond series of values for the conditional variable
	 * @return series of local conditional MI values
	 */
	public double[] computeLocalFromPreviousObservations(int var1[], int var2[], int cond[]){
		int rows = var1.length;

		double[] localCondMi = new double[rows];
		average = 0;
		max = 0;
		min = 0;

		int var1Val, var2Val, condVal;
		double logTerm;
		for (int r = 0; r < rows; r++) {
			var1Val = var1[r];
			var2Val = var2[r];
			condVal = cond[r];
			// Now compute the local value
			logTerm = ((double) firstSecondCondCount[var1Val][var2Val][condVal] / (double) firstCondCount[var1Val][condVal]) /
		 		((double) secondCondCount[var2Val][condVal] / (double) condCount[condVal]);
			localCondMi[r] = Math.log(logTerm) / log_2;
			average += localCondMi[r];
			if (localCondMi[r] > max) {
				max = localCondMi[r];
			} else if (localCondMi[r] < min) {
				min = localCondMi[r];
			}
		}

		average = average/(double) rows;
		condMiComputed = true;

		return localCondMi;
	}

	/**
	 * Generate a bootstrapped distribution of what the conditional MI would look like,
	 * under a null hypothesis that the source values of our
	 * samples had no relation to the destination value (in the
	 * context of the conditional).
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for MI,
	 * conditional MI and TE.
	 * <b>Note that this method currently fixes the relationship
	 * between variable 2 and the conditional, and shuffles
	 * variable 1 with respect to these.</b>
	 * </p>
	 * 
	 * <p>Note that if several disjoint time-series have been added 
	 * as observations using {@link #addObservations(int[], int[])} etc.,
	 * then these separate "trials" will be mixed up in the generation
	 * of surrogates here.</p>
	 * 
	 * @param numPermutationsToCheck number of surrogate samples to bootstrap
	 *  to generate the distribution.
	 * @return the distribution of conditional MI scores under this null hypothesis.
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
	 * Generate a bootstrapped distribution of what the conditional MI would look like,
	 * under a null hypothesis that the source values of our
	 * samples had no relation to the destination values.
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for  
	 * a conditional mutual information. Basically, the marginal PDFs
	 * of each marginal
	 * are preserved, while their joint PDF is destroyed, and the 
	 * distribution of conditional MI under these conditions is generated.
	 * <b>Note that this method currently fixes the relationship
	 * between variable 2 and the conditional, and shuffles
	 * variable 1 with respect to these.</b>
	 * </p>
	 * TODO Need to alter the method signature to allow callers to specify
	 * which variable is shuffled. (Note to self: when doing this, will
	 * need to update machine learning code to the new method signature)
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
	 * @param newOrderings a specification of how to shuffle the values
	 *  of variable 1
	 *  to create the surrogates to generate the distribution with. The first
	 *  index is the permutation number (i.e. newOrderings.length is the number
	 *  of surrogate samples we use to bootstrap to generate the distribution here.)
	 *  Each array newOrderings[i] should be an array of length N (where
	 *  would be the value returned by {@link #getNumObservations()}),
	 *  containing a permutation of the values in 0..(N-1).
	 * @return the distribution of conditional MI scores under this null hypothesis.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception where the length of each permutation in newOrderings
	 *   is not equal to the number N samples that were previously supplied.
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) {
		
		double actualCondMI = computeAverageLocalOfObservations();
		
		int numPermutationsToCheck = newOrderings.length;
		
		// Reconstruct the observed values of the variables in some order
		int[] var1Values = new int[observations];
		int[] var2Values = new int[observations];
		int[] condValues = new int[observations];
		int t_s = 0;
		for (int val1 = 0; val1 < base1; val1++) {
			for (int val2 = 0; val2 < base2; val2++) {
				for (int condVal = 0; condVal < condBase; condVal++) {
					int numberOfSamples = firstSecondCondCount[val1][val2][condVal];
					MatrixUtils.fill(var1Values, val1, t_s, numberOfSamples);
					MatrixUtils.fill(var2Values, val2, t_s, numberOfSamples);
					MatrixUtils.fill(condValues, condVal, t_s, numberOfSamples);
					t_s += numberOfSamples;
				}
			}
		}
		// We now have arrays of the values that were observed for each
		//  variable, in a random order (well, actually, in order of
		//  increasing joint value of the observations, but this doesn't
		//  matter because:). We will now extract randomly ordered
		//  time series of var1Values, to bootstrap the distribution
		//  of conditional MI values under the null hypothesis.
				
		ConditionalMutualInformationCalculatorDiscrete condMi2 =
				new ConditionalMutualInformationCalculatorDiscrete(base1, base2, condBase);
		condMi2.initialise();
		// Set up the joint counts which remain the same under reordering
		//  of variable 1:
		condMi2.observations = observations;
		condMi2.secondCondCount = secondCondCount;
		condMi2.condCount = condCount;
		int countWhereMIIsMoreSignificantThanOriginal = 0;
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		for (int p = 0; p < numPermutationsToCheck; p++) {
			// Generate a new re-ordered data set for the 1st variable
			int[] newData1 = MatrixUtils.extractSelectedTimePoints(var1Values, newOrderings[p]);
			// Compute the required joint probability distributions:
			MatrixUtils.fill(condMi2.firstCondCount, 0);
			MatrixUtils.fill(condMi2.firstSecondCondCount, 0);
			for (int t = 0; t < observations; t++) {
				condMi2.firstCondCount[newData1[t]][condValues[t]]++;
				condMi2.firstSecondCondCount[newData1[t]][var2Values[t]][condValues[t]]++;
			}
			// And get a cond MI value for this realisation of var1Values:
			double newCondMI = condMi2.computeAverageLocalOfObservations();
			measDistribution.distribution[p] = newCondMI;
			if (newCondMI >= actualCondMI) {
				countWhereMIIsMoreSignificantThanOriginal++;
			}
		}
		
		// And return the significance
		measDistribution.pValue = (double) countWhereMIIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = actualCondMI;
		return measDistribution;
	}

	@Override
	public AnalyticMeasurementDistribution computeSignificance() {
		if (!condMiComputed) {
			computeAverageLocalOfObservations();
		}
		return new ChiSquareMeasurementDistribution(2.0*((double)observations)*average,
				(base1 - 1)*(base2 - 1)*condBase);
	}

	/**
	 * Standalone routine to 
	 * compute local conditional MI for given variables.
	 * Return a temporal array of local values.
	 * 
	 * @param var1 series of values for the first variable
	 * @param var2 series of values for the second variable
	 * @param cond series of values for the conditional variable
	 * @return series of local conditional MI values
	 */
	public double[] computeLocal(int var1[], int var2[], int cond[]) {
		
		initialise();
		addObservations(var1, var2, cond);
		return computeLocalFromPreviousObservations(var1, var2, cond);
	}

	/**
	 * Standalone routine to 
	 * compute average conditional MI for a given set of variables.
	 * 
	 * @param var1 series of values for the first variable
	 * @param var2 series of values for the second variable
	 * @param cond series of values for the conditional variable
	 * @return average conditional MI for these values
	 */
	public double computeAverageLocal(int var1[], int var2[], int cond[]) {
		
		initialise();
		addObservations(var1, var2, cond);
		return computeAverageLocalOfObservations();
	}
}
