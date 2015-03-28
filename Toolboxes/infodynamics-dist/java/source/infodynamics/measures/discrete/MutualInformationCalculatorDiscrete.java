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
 * <p>Mutual information (MI) calculator for univariate discrete (int[]) data.</p>
 * 
 * <p>Usage of the class is intended to follow this paradigm:</p>
 * <ol>
 * 		<li>Construct the calculator: {@link #MutualInformationCalculatorDiscrete(int)}
 * 			or {@link #MutualInformationCalculatorDiscrete(int, int)};</li>
 *		<li>Initialise the calculator using {@link #initialise()};</li>
 * 		<li>Provide the observations/samples for the calculator
 *      	to set up the PDFs, using one or more calls to
 * 			sets of {@link #addObservations(int[], int[])} methods, then</li>
 * 		<li>Compute the required quantities, being one or more of:
 * 			<ul>
 * 				<li>the average MI: {@link #computeAverageLocalOfObservations()};</li>
 * 				<li>local MI values, such as
 * 				{@link #computeLocalFromPreviousObservations(int[], int[])};</li>
 * 				<li>comparison to null distribution, such as
 * 				{@link #computeSignificance()};</li>
 * 				<li>and variants of these.</li>
 * 			</ul>
 * 		</li>
 * 		<li>As an alternative to steps 3 and 4, the user may undertake
 * 			standalone computation from a single set of observations, via
 *  		e.g.: {@link #computeLocal(int[][], int, int)}.</li>
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
public class MutualInformationCalculatorDiscrete extends InfoMeasureCalculatorDiscrete 
	implements ChannelCalculatorDiscrete, AnalyticNullDistributionComputer {

	private int timeDiff = 0;
	private int[][]	jointCount = null; // Count for (i[t-timeDiff], j[t]) tuples
	private int[] iCount = null; // Count for i[t-timeDiff]		
	private int[] jCount = null; // Count for j[t]

	protected boolean miComputed = false;
	
	/**
	 * Construct a new MI calculator with default time difference of 0
	 *  between the variables
	 * 
	 * @param base number of symbols for each variable.
	 *        E.g. binary variables are in base-2.
	 * @throws Exception
	 */
	public MutualInformationCalculatorDiscrete(int base) throws Exception {
		this(base, 0);
	}
	
	/**
	 * Create a new mutual information calculator
	 * 
	 * @param base number of symbols for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param timeDiff number of time steps across which to compute
	 *   MI for given time series
	 * @throws Exception when timeDiff < 0
	 */
	public MutualInformationCalculatorDiscrete(int base, int timeDiff) throws Exception {
		super(base);
		if (timeDiff < 0) {
			throw new Exception("timeDiff must be >= 0");
		}
		this.timeDiff = timeDiff;
		jointCount = new int[base][base];
		iCount = new int[base];
		jCount = new int[base];
	}

	@Override
	public void initialise(){
		super.initialise();
		miComputed = false;		
		MatrixUtils.fill(iCount, 0);
		MatrixUtils.fill(jCount, 0);
		MatrixUtils.fill(jointCount, 0);
	}
	
	/**
	 * {@inheritDoc}
	 * 
	 * Pairs for MI are between the arrays var1 and var2, separated in time by timeDiff
	 * (var1 is first).
	 *
	 */
	@Override
	public void addObservations(int[] var1, int[] var2) {
		int timeSteps = var1.length;
		// int columns = states[0].length;
		// increment the count of observations:
		observations += (timeSteps - timeDiff); 
		
		// 1. Count the tuples observed
		int iVal, jVal;
		for (int r = timeDiff; r < timeSteps; r++) {
			// Add to the count for this particular pair:
			iVal = var1[r-timeDiff];
			jVal = var2[r];
			jointCount[iVal][jVal]++;
			iCount[iVal]++;
			jCount[jVal]++;
		}
	}

	/**
	 * {@inheritDoc}
	 * 
	 * Pairs for MI are between columns iCol and jCol, separated in time by timeDiff (i is first).
	 *
	 */
	@Override
	public void addObservations(int states[][], int iCol, int jCol) {
		int rows = states.length;
		// int columns = states[0].length;
		// increment the count of observations:
		observations += (rows - timeDiff); 
		
		// 1. Count the tuples observed
		int iVal, jVal;
		for (int r = timeDiff; r < rows; r++) {
			// Add to the count for this particular pair:
			iVal = states[r-timeDiff][iCol];
			jVal = states[r][jCol];
			jointCount[iVal][jVal]++;
			iCount[iVal]++;
			jCount[jVal]++;
		}
	}
	
	@Override
	public double computeAverageLocalOfObservations() {
		double mi = 0.0;
		double miCont = 0.0;

		max = 0;
		min = 0;
		double meanSqLocals = 0;
		if (debug) {
			System.out.println("i\tj\tp_i\tp_j\tp_joint\tlocal");
		}
		for (int i = 0; i < base; i++) {
			// compute p_i
			double probi = (double) iCount[i] / (double) observations;
			for (int j = 0; j < base; j++) {
				// compute p_j
				double probj = (double) jCount[j] / (double) observations;
				// compute p(veci=i, vecj=j)
				double jointProb = (double) jointCount[i][j] / (double) observations;
				// Compute MI contribution:
				if (jointProb * probi * probj > 0.0) {
					double localValue = Math.log(jointProb / (probi * probj)) / log_2;
					miCont = jointProb * localValue;
					if (debug) {
						System.out.printf("%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n",
								i, j, probi, probj, jointProb, localValue);
					}
					if (localValue > max) {
						max = localValue;
					} else if (localValue < min) {
						min = localValue;
					}
					// Add this contribution to the mean 
					//  of the squared local values
					meanSqLocals += miCont * localValue;
				} else {
					miCont = 0.0;
				}
				mi += miCont;
			}
		}
		
		average = mi;
		miComputed = true;
		std = Math.sqrt(meanSqLocals - average * average);

		return mi;
	}

	@Override
	public EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) {
		RandomGenerator rg = new RandomGenerator();
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(observations, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}
	
	/**
	 * Generate a bootstrapped distribution of what the MI would look like,
	 * under a null hypothesis that the source values of our
	 * samples had no relation to the destination values.
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for  
	 * a mutual information. Basically, the marginal PDFs
	 * of each marginal
	 * are preserved, while their joint PDF is destroyed, and the 
	 * distribution of MI under these conditions is generated.</p>
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
	 * @return the distribution of MI scores under this null hypothesis.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception where the length of each permutation in newOrderings
	 *   is not equal to the number N samples that were previously supplied.
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) {
		double actualMI = computeAverageLocalOfObservations();
		
		int numPermutationsToCheck = newOrderings.length;
		
		// Reconstruct the values of the first and second variables (not necessarily in order)
		int[] iValues = new int[observations];
		int[] jValues = new int[observations];
		int t_i = 0;
		int t_j = 0;
		for (int iVal = 0; iVal < base; iVal++) {
			int numberOfSamplesI = iCount[iVal];
			MatrixUtils.fill(iValues, iVal, t_i, numberOfSamplesI);
			t_i += numberOfSamplesI;
			int numberOfSamplesJ = jCount[iVal];
			MatrixUtils.fill(jValues, iVal, t_j, numberOfSamplesJ);
			t_j += numberOfSamplesJ;
		}
		
		MutualInformationCalculatorDiscrete mi2;
		try {
			mi2 = new MutualInformationCalculatorDiscrete(base, timeDiff);
		} catch (Exception e) {
			// The only possible exception is if timeDiff < 0, which 
			// it cannot be. Shut down the JVM
			throw new Error("timeDiff parameter took on value < 0 after being checked at construction");
		}
		mi2.initialise();
		mi2.observations = observations;
		mi2.iCount = iCount;
		mi2.jCount = jCount;
		int countWhereMIIsMoreSignificantThanOriginal = 0;
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		for (int p = 0; p < numPermutationsToCheck; p++) {
			// Generate a new re-ordered data set for the i variable
			int[] newDataI = MatrixUtils.extractSelectedTimePoints(iValues, newOrderings[p]);
			// compute the joint probability distribution
			MatrixUtils.fill(mi2.jointCount, 0);
			for (int t = 0; t < observations; t++) {
				mi2.jointCount[newDataI[t]][jValues[t]]++;
			}
			// And get an MI value for this realisation:
			double newMI = mi2.computeAverageLocalOfObservations();
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

	@Override
	public AnalyticMeasurementDistribution computeSignificance() {
		if (!miComputed) {
			computeAverageLocalOfObservations();
		}
		return new ChiSquareMeasurementDistribution(2.0*((double)observations)*average,
				(base - 1) * (base - 1));
	}
	
	/**
	 * Computes local mutual information (or pointwise mutual information)
	 *  for the given (single) specific values, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 *  
	 * @param val1 single specific value of variable 1
	 * @param val2 single specific value of variable 2
	 * @return a local mutual information value for this pair of observations
	 * @throws Exception if this pair were not observed together in the
	 *  previously supplied observations
	 */
	public double computeLocalFromPreviousObservations(int val1, int val2) throws Exception{
		
		double logTerm = ( (double) jointCount[val1][val2] ) /
			  		  ( (double) jCount[val2] *
			  			(double) iCount[val1] );
		// Now account for the fact that we've
		//  just used counts rather than probabilities,
		//  and we've got two counts on the bottom
		//  but one count on the top:
		logTerm *= (double) observations;
		double localMI = Math.log(logTerm) / log_2;
		
		return localMI;
	}
	
	/**
	 * Computes local mutual information (or pointwise mutual information)
	 *  for the given states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 *  
	 * @param var1 new states of variable 1
	 * @param var2 new states of variable 2. Should be same length as var1
	 * @return array of local mutual information values for each
	 *  observation of (var1, var2). Note - if timeDiff > 0, then the
	 *  return length will be var1.length - timeDiff. 
	 */
	public double[] computeLocalFromPreviousObservations(int[] var1, int[] var2) throws Exception{
		
		if (var1.length != var2.length) {
			throw new Exception("var1 and var2 must have the same number of observations");
		}
		double[] localMI = new double[var1.length - timeDiff];
		
		double logTerm = 0.0;
		for (int r = timeDiff; r < var1.length; r++) {
			int iVal = var1[r-timeDiff];
			int jVal = var2[r];
			logTerm = ( (double) jointCount[iVal][jVal] ) /
			  		  ( (double) jCount[jVal] *
			  			(double) iCount[iVal] );
			// Now account for the fact that we've
			//  just used counts rather than probabilities,
			//  and we've got two counts on the bottom
			//  but one count on the top:
			logTerm *= (double) observations;
			localMI[r] = Math.log(logTerm) / log_2;
			average += localMI[r];
			if (localMI[r] > max) {
				max = localMI[r];
			} else if (localMI[r] < min) {
				min = localMI[r];
			}
		}
		average = average/(double) observations;
		miComputed = true;

		return localMI;
	}
	
	/**
	 * Computes local mutual information (or pointwise mutual information)
	 *  for the given states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 *  
	 * @param states 2D time series of observations (first index time,
	 *  second is variable index)
	 * @param iCol column number for first variable
	 * @param jCol column number for second variable
	 * @return array of local mutual information values for each
	 *  observation of (var1, var2). Note - if timeDiff > 0, then the
	 *  return length will be var1.length - timeDiff. 
	 */
	public double[] computeLocalFromPreviousObservations(int states[][], int iCol, int jCol){
		int rows = states.length;
		//int columns = states[0].length;

		// Allocate for all rows even though we'll leave the first ones as zeros
		double[] localMI = new double[rows];
		int iVal, jVal;
		double logTerm = 0.0;
		for (int r = timeDiff; r < rows; r++) {
			iVal = states[r-timeDiff][iCol];
			jVal = states[r][jCol];
			logTerm = ( (double) jointCount[iVal][jVal] ) /
			  		  ( (double) jCount[jVal] *
			  			(double) iCount[iVal] );
			// Now account for the fact that we've
			//  just used counts rather than probabilities,
			//  and we've got two counts on the bottom
			//  but one count on the top:
			logTerm *= (double) observations;
			localMI[r] = Math.log(logTerm) / log_2;
			average += localMI[r];
			if (localMI[r] > max) {
				max = localMI[r];
			} else if (localMI[r] < min) {
				min = localMI[r];
			}
		}
		average = average/(double) observations;
		
		return localMI;
		
	}
	
	/**
	 * Standalone routine to 
	 * compute local mutual information (or pointwise mutual information)
	 *  across a 2D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros
	 * 
	 * @param states 2D time series of observations (first index time,
	 *  second is variable index)
	 * @param iCol column number for first variable
	 * @param jCol column number for second variable
	 * @return array of local mutual information values for each
	 *  observation of (var1, var2). Note - if timeDiff > 0, then the
	 *  return length will be var1.length - timeDiff. 
	 */
	public double[] computeLocal(int states[][], int iCol, int jCol) {
		initialise();
		addObservations(states, iCol, jCol);
		return computeLocalFromPreviousObservations(states, iCol, jCol);
		
	}
}
