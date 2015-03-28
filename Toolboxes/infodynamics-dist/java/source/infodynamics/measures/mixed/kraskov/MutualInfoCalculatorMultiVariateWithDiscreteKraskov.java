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

package infodynamics.measures.mixed.kraskov;

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.measures.mixed.MutualInfoCalculatorMultiVariateWithDiscrete;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

/**
 * <p>Compute the Mutual Information between a vector of continuous variables and discrete
 *  variable using the Kraskov estimation method.</p>
 * <p>Uses Kraskov method type 2, since type 1 only looks at points with
 * distances strictly less than the kth variable, which won't work for one marginal
 * being discrete.</p>
 * <p>I have noticed that there are quite large bias negative values here
 * where small K is used (e.g. for binary data splitting continuous into
 * two distinct groups, 1400 observations, K=4 has bias ~ -.35)</p>
 * 
 * <p>These calculators are <b>EXPERIMENTAL</b> -- not properly tested,
 * and not well documented. The intended calling pattern is similar to
 * {@link MutualInfoCalculatorMultiVariate}
 * </p>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class MutualInfoCalculatorMultiVariateWithDiscreteKraskov implements MutualInfoCalculatorMultiVariateWithDiscrete {

	// Multiplier used in hueristic for determining whether to use a linear search
	//  for min kth element or a binary search.
	protected static final double CUTOFF_MULTIPLIER = 1.5;
	
	/**
	 * we compute distances to the kth neighbour
	 */
	protected int k = 4;
	protected double[][] continuousData;
	protected int[] discreteData;
	protected int[] counts;
	/**
	 * Number of possible states of the discrete variable
	 */
	protected int base;
	/**
	 * Number of dimenions of the joint continuous variable
	 */
	protected int dimensions;
	protected boolean debug;
	protected double mi;
	protected boolean miComputed;
	
	protected EuclideanUtils normCalculator;
	// Storage for the norms from each observation to each other one
	protected double[][] xNorms;
	// Keep the norms each time (making reordering very quick)
	//  (Should only be set to false for testing)
	public static boolean tryKeepAllPairsNorms = true;
	public static int MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM = 2000;
	
	public final static String PROP_K = "k";
	public final static String PROP_NORM_TYPE = "NORM_TYPE";	
	public static final String PROP_NORMALISE = "NORMALISE";
	
	/**
	 * Track whether we're going to normalise the joint variables individually
	 */
	protected boolean normalise = true;
	/**
	 * Track the means of the joint variables if we are normalising them
	 */
	protected double[] means;
	/**
	 * Track the std devs of the joint variables if we are normalising them
	 */
	protected double[] stds;

	public MutualInfoCalculatorMultiVariateWithDiscreteKraskov() {
		super();
		normCalculator = new EuclideanUtils(EuclideanUtils.NORM_MAX_NORM);
	}

	/**
	 * Initialise the calculator.
	 * 
	 * @param dimensions number of joint continuous variables
	 * @param base number of discrete states
	 */
	public void initialise(int dimensions, int base) {
		mi = 0.0;
		miComputed = false;
		xNorms = null;
		continuousData = null;
		means = null;
		stds = null;
		discreteData = null;
		this.dimensions = dimensions;
		this.base = base;
	}

	/**
	 * Sets properties for the calculator.
	 * Valid properties include:
	 * <ul>
	 *  <li>{@link #PROP_K} - number of neighbouring points in joint kernel space (default 4)</li>
	 * 	<li>{@link #PROP_NORM_TYPE}</li> - normalization type to apply to 
	 * 		working out the norms between the points in each marginal space.
	 * 		Options are defined by {@link EuclideanUtils#setNormToUse(String)} -
	 * 		default is {@link EuclideanUtils#NORM_MAX_NORM}.
	 *  <li>{@link #PROP_NORMALISE} - whether to normalise the individual
	 *      variables (true by default)</li>
	 * </ul>
	 * 
	 * @param propertyName name of the property to set
	 * @param propertyValue value to set on that property
	 */
	public void setProperty(String propertyName, String propertyValue) {
		if (propertyName.equalsIgnoreCase(PROP_K)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
			normCalculator.setNormToUse(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORMALISE)) {
			normalise = Boolean.parseBoolean(propertyValue);
		}
	}

	public void addObservations(double[][] source, double[][] destination) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void addObservations(double[][] source, double[][] destination, int startTime, int numTimeSteps) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] source, double[][] destination, boolean[] sourceValid, boolean[] destValid) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] source, double[][] destination, boolean[][] sourceValid, boolean[][] destValid) throws Exception {
		throw new RuntimeException("Not implemented yet");
	}

	public void startAddObservations() {
		throw new RuntimeException("Not implemented yet");
	}

	public void finaliseAddObservations() {
		throw new RuntimeException("Not implemented yet");
	}

	public void setObservations(double[][] continuousObservations,
			int[] discreteObservations) throws Exception {
		if (continuousObservations.length != discreteObservations.length) {
			throw new Exception("Time steps for observations2 " +
					discreteObservations.length + " does not match the length " +
					"of observations1 " + continuousObservations.length);
		}
		if (continuousObservations[0].length == 0) {
			throw new Exception("Computing MI with a null set of data");
		}
		if (continuousObservations[0].length != dimensions) {
			throw new Exception("The continuous observations do not have the expected number of variables (" + dimensions + ")");
		}
		continuousData = continuousObservations;
		discreteData = discreteObservations;
		if (normalise) {
			// Take a copy since we're going to normalise it
			// And we need to keep the means/stds ready to normalise local values
			//  that are supplied later:
			means = MatrixUtils.means(continuousObservations);
			stds = MatrixUtils.stdDevs(continuousObservations, means);
			continuousData = MatrixUtils.normaliseIntoNewArray(continuousObservations, means, stds);
		}
		// count the discrete states:
		counts = new int[base];
		for (int t = 0; t < discreteData.length; t++) {
			counts[discreteData[t]]++;
		}
		for (int b = 0; b < counts.length; b++) {
			if (counts[b] < k) {
				throw new RuntimeException("This implementation assumes there are at least k items in each discrete bin");
			}
		}
	}

	/**
	 * Compute the norms for each marginal time series
	 *
	 */
	protected void computeNorms() {
		int N = continuousData.length; // number of observations
		
		xNorms = new double[N][N];
		for (int t = 0; t < N; t++) {
			// Compute the norms from t to all other time points
			for (int t2 = 0; t2 < N; t2++) {
				if (t2 == t) {
					xNorms[t][t2] = Double.POSITIVE_INFINITY;
					continue;
				}
				// Compute norm in the continuous space
				xNorms[t][t2] = normCalculator.norm(continuousData[t], continuousData[t2]);
			}
		}
	}
	
	/**
	 * Compute what the average MI would look like were the second time series reordered
	 *  as per the array of time indices in reordering.
	 * The user should ensure that all values 0..N-1 are represented exactly once in the
	 *  array reordering and that no other values are included here. 
	 * 
	 * @param reordering
	 * @return
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations(int[] reordering) throws Exception {
		int N = continuousData.length; // number of observations
		if (!tryKeepAllPairsNorms || (N > MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM)) {
			// Generate a new re-ordered set of discrete data
			int[] originalDiscreteData = discreteData;
			discreteData = MatrixUtils.extractSelectedTimePoints(discreteData, reordering);
			// Compute the MI
			double newMI = computeAverageLocalOfObservationsWhileComputingDistances();
			// restore data2
			discreteData = originalDiscreteData;
			return newMI;
		}
		
		// Otherwise we will use the norms we've already computed, and use a "virtual"
		//  reordered data2.
		int[] reorderedDiscreteData = MatrixUtils.extractSelectedTimePoints(discreteData, reordering);

		if (xNorms == null) {
			computeNorms();
		}

		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double avNx = 0;
		double avNy = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps_x and eps_y for this time step:
			//  using x norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			// Then find the k closest neighbours in the same discrete bin
			double eps_x = MatrixUtils.kthMinSubjectTo(xNorms[t], k, reorderedDiscreteData, reorderedDiscreteData[t]);			
			
			// Count the number of points whose x distance is less
			//  than or equal to eps_x
			int n_x = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (xNorms[t][t2] <= eps_x) {
					n_x++;
				}
			}
			// n_y is number of points in that bin minus that point
			int n_y = counts[reorderedDiscreteData[t]] - 1;
			avNx += n_x;
			avNy += n_y;
			// And take the digamma before adding into the 
			//  average:
			averageDiGammas += MathsUtils.digamma(n_x) + MathsUtils.digamma(n_y);
		}
		averageDiGammas /= (double) N;
		if (debug) {
			avNx /= (double)N;
			avNy /= (double)N;
			System.out.println(String.format("Average n_x=%.3f, Average n_y=%.3f", avNx, avNy));
		}
		
		mi = MathsUtils.digamma(k) - 1.0/(double)k - averageDiGammas + MathsUtils.digamma(N);
		miComputed = true;
		return mi;
	}

	public double computeAverageLocalOfObservations() throws Exception {
		if (!tryKeepAllPairsNorms || (continuousData.length > MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM)) {
			return computeAverageLocalOfObservationsWhileComputingDistances();
		}
		
		// Postcondition: we'll compute the norms before the main loop,
		//  unless they have already been computed:
		
		if (xNorms == null) {
			computeNorms();
		}
		int N = continuousData.length; // number of observations

		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double avNx = 0;
		double avNy = 0;
		
		double testSum = 0.0; // Used for debugging prints
		
		for (int t = 0; t < N; t++) {
			// Compute eps_x and eps_y for this time step:
			//  using x norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			// Then find the k closest neighbours in the same discrete bin
			double eps_x = MatrixUtils.kthMinSubjectTo(xNorms[t], k, discreteData, discreteData[t]);			

			// Count the number of points whose x distance is less
			//  than or equal to eps_x (not including this point)
			int n_x = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if ((t2 != t) && (xNorms[t][t2] <= eps_x)) {
					// xNorms has infinity for t == t2 anyway ...
					n_x++;
				}
			}
			// n_y is number of points in that bin, minus this point
			int n_y = counts[discreteData[t]] - 1;
			avNx += n_x;
			avNy += n_y;
			// And take the digamma before adding into the 
			//  average:
			double localSum = MathsUtils.digamma(n_x) + MathsUtils.digamma(n_y);
			averageDiGammas += localSum;

			if (debug) {
				double localValue = MathsUtils.digamma(k) - 1.0/(double)k - localSum + MathsUtils.digamma(N);
				testSum += localValue;
				if (dimensions == 1) {
					System.out.printf("t=%d: x=%.3f, eps_x=%.3f, n_x=%d, n_y=%d, local=%.3f, running total = %.5f\n",
							t, continuousData[t][0], eps_x, n_x, n_y, localValue, testSum);
				} else {
					System.out.printf("t=%d: eps_x=%.3f, n_x=%d, n_y=%d, local=%.3f, running total = %.5f\n",
							t, eps_x, n_x, n_y, localValue, testSum);
				}
			}
		}
		averageDiGammas /= (double) N;
		if (debug) {
			avNx /= (double)N;
			avNy /= (double)N;
			System.out.println(String.format("Average n_x=%.3f (-> digam=%.3f %.3f), Average n_y=%.3f (-> digam=%.3f)",
					avNx, MathsUtils.digamma((int) avNx), MathsUtils.digamma((int) avNx - 1), avNy, MathsUtils.digamma((int) avNy)));
			System.out.printf("Independent average num in joint box is %.3f\n", (avNx * avNy / (double) N));
			System.out.println(String.format("digamma(k)=%.3f - 1/k=%.3f - averageDiGammas=%.3f + digamma(N)=%.3f\n",
					MathsUtils.digamma(k), 1.0/(double)k, averageDiGammas, MathsUtils.digamma(N)));
		}
		
		mi = MathsUtils.digamma(k) - 1.0/(double)k - averageDiGammas + MathsUtils.digamma(N);
		miComputed = true;
		return mi;
	}

	/**
	 * This method correctly computes the average local MI, but recomputes the x and y 
	 *  distances between all tuples in time.
	 * Kept here for cases where we have too many observations
	 *  to keep the norm between all pairs, and for testing purposes.
	 * 
	 * @return
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservationsWhileComputingDistances() throws Exception {
		int N = continuousData.length; // number of observations

		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double avNx = 0;
		double avNy = 0;
		
		double testSum = 0.0; // Used for debugging prints
		
		for (int t = 0; t < N; t++) {
			// Compute eps_x and eps_y for this time step:
			//  First get x norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			double[] norms = new double[N]; 
			for (int t2 = 0; t2 < N; t2++) {
				if (t2 == t) {
					norms[t2] = Double.POSITIVE_INFINITY;
					continue;
				}
				// Compute norm in the continuous space
				norms[t2] = normCalculator.norm(continuousData[t], continuousData[t2]);
			}

			// Then find the k closest neighbours in the same discrete bin
			double eps_x = MatrixUtils.kthMinSubjectTo(norms, k, discreteData, discreteData[t]);			

			// Count the number of points whose x distance is less
			//  than or equal to eps_x
			int n_x = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (norms[t2] <= eps_x) {
					n_x++;
				}
			}
			// n_y is number of points in that bin, minus that point
			int n_y = counts[discreteData[t]] - 1;
			avNx += n_x;
			avNy += n_y;
			// And take the digamma before adding into the 
			//  average:
			double localSum = MathsUtils.digamma(n_x) + MathsUtils.digamma(n_y);
			averageDiGammas += localSum;

			if (debug) {
				double localValue = MathsUtils.digamma(k) - 1.0/(double)k - localSum + MathsUtils.digamma(N);
				testSum += localValue;
				if (dimensions == 1) {
					System.out.printf("t=%d: x=%.3f, eps_x=%.3f, n_x=%d, n_y=%d, local=%.3f, running total = %.5f\n",
							t, continuousData[t][0], eps_x, n_x, n_y, localValue, testSum);
				} else {
					System.out.printf("t=%d: eps_x=%.3f, n_x=%d, n_y=%d, local=%.3f, running total = %.5f\n",
							t, eps_x, n_x, n_y, localValue, testSum);
				}
			}
		}
		averageDiGammas /= (double) N;
		if (debug) {
			avNx /= (double)N;
			avNy /= (double)N;
			System.out.println(String.format("Average n_x=%.3f, Average n_y=%.3f", avNx, avNy));
		}
		
		mi = MathsUtils.digamma(k) - 1.0/(double)k - averageDiGammas + MathsUtils.digamma(N);
		miComputed = true;
		return mi;
	}

	/**
	 * Compute the significance of the mutual information of the previously supplied observations.
	 * We destroy the p(x,y) correlations, while retaining the p(x), p(y) marginals, to check how
	 *  significant this mutual information actually was.
	 *  
	 * This is in the spirit of Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128
	 *  which was performed for Transfer entropy.
	 * 
	 * @param numPermutationsToCheck
	 * @return the proportion of MI scores from the distribution which have higher or equal MIs to ours.
	 *  (i.e. 1 - CDF of our score)
	 */
	public synchronized EmpiricalMeasurementDistribution computeSignificance(int numPermutationsToCheck) throws Exception {
		// Generate the re-ordered indices:
		RandomGenerator rg = new RandomGenerator();
		// (Not necessary to check for distinct random perturbations)
		int[][] newOrderings = rg.generateRandomPerturbations(continuousData.length, numPermutationsToCheck);
		return computeSignificance(newOrderings);
	}
	
	/**
	 * Compute the significance of the mutual information of the previously supplied observations.
	 * We destroy the p(x,y) correlations, while retaining the p(x), p(y) marginals, to check how
	 *  significant this mutual information actually was.
	 *  
	 * This is in the spirit of Chavez et. al., "Statistical assessment of nonlinear causality:
	 *  application to epileptic EEG signals", Journal of Neuroscience Methods 124 (2003) 113-128
	 *  which was performed for Transfer entropy.
	 * 
	 * @param newOrderings the specific new orderings to use
	 * @return the proportion of MI scores from the distribution which have higher or equal MIs to ours.
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings) throws Exception {
		
		int numPermutationsToCheck = newOrderings.length;
		if (!miComputed) {
			computeAverageLocalOfObservations();
		}
		// Store the real observations and their MI:
		double actualMI = mi;
		
		EmpiricalMeasurementDistribution measDistribution = new EmpiricalMeasurementDistribution(numPermutationsToCheck);
		
		int countWhereMiIsMoreSignificantThanOriginal = 0;
		for (int i = 0; i < numPermutationsToCheck; i++) {
			// Compute the MI under this reordering
			double newMI = computeAverageLocalOfObservations(newOrderings[i]);
			measDistribution.distribution[i] = newMI;
			if (debug){
				System.out.println("New MI was " + newMI);
			}
			if (newMI >= actualMI) {
				countWhereMiIsMoreSignificantThanOriginal++;
			}
		}
		
		// Restore the actual MI and the observations
		mi = actualMI;

		// And return the significance
		measDistribution.pValue = (double) countWhereMiIsMoreSignificantThanOriginal / (double) numPermutationsToCheck;
		measDistribution.actualValue = mi;
		return measDistribution;
	}

	/**
	 * <p>Compute the local MI values for the given observations,
	 *  using the previously set observations to compute the PDFs.
	 *  That is to say, we will evaluate the required counts for each
	 *  observation here based on the state space constructed from the
	 *  previous observations only (not these ones). This is non-standard,
	 *  and is not considered in the Kraskov paper. It is slightly unclear
	 *  how to account for the fact that this observation itself is not
	 *  in the data set used for computing the PDFs - I think though that
	 *  it should be ignored, since the data point was not counted in its
	 *  own counts in the standard version anyway.</p>
	 *  
	 *  <p><b>Importantly</b>, the supplied observations are intended to be new observations,
	 *  not those fed in to compute the PDFs from. There would be
	 *  some subtle changes necessary to accomodate computing locals on
	 *  the same data set (e.g. not counting the current point as one
	 *  of those within eps_x etc.). 
	 *  </p>
	 *  
	 *  @param continuousStates multivariate observations of the continuous variable
	 *    (1st index is time, 2nd is variable number)
	 *  @param discreteStates unvariate observations of the discrete variable
	 * 
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] continuousNewStates,
			int[] discreteNewStates) throws Exception {
		
		if (normalise) {
			// The stored observations continuousData have been normalised
			//  according to their stored means and stds; we need to 
			//  normalise the incoming observations the same way before
			//  comparing them
			continuousNewStates = MatrixUtils.normaliseIntoNewArray(
					continuousNewStates, means, stds);
		}
		
		int N = continuousNewStates.length; // number of observations
		double[] locals = new double[N];
		
		double fixedPartOfLocals = MathsUtils.digamma(k) - 1.0/(double)k +
									MathsUtils.digamma(N);
		double testSum = 0.0;
		if (debug) {
			System.out.printf("digamma(k)=%.3f - 1/k=%.3f + digamma(N)=%.3f\n",
					MathsUtils.digamma(k), 1.0/(double)k, MathsUtils.digamma(N));
		}
		double avNx = 0;
		double avNy = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps_x and eps_y for this time step:
			//  First get x norms to all points in the previously
			//   given observations.
			double[] norms = new double[continuousData.length];
			for (int t2 = 0; t2 < continuousData.length; t2++) {
				// Compute norm in the continuous space
				norms[t2] = normCalculator.norm(continuousNewStates[t], continuousData[t2]);
			}

			// Then find the k closest neighbours in the same discrete bin
			double eps_x = MatrixUtils.kthMinSubjectTo(norms, k, discreteData, discreteNewStates[t]);			

			// Count the number of points whose x distance is less
			//  than or equal to eps_x
			int n_x = 0;
			for (int t2 = 0; t2 < continuousData.length; t2++) {
				if (norms[t2] <= eps_x) {
					n_x++;
				}
			}
			// n_y is number of points in that discrete bin
			int n_y = counts[discreteData[t]];
			avNx += n_x;
			avNy += n_y;
			// Now compute the local value:
			locals[t] = fixedPartOfLocals -
					MathsUtils.digamma(n_x) - MathsUtils.digamma(n_y);
			if (debug) {
				testSum += locals[t];
				if (dimensions == 1) {
					System.out.printf("t=%d: x=%.3f, eps_x=%.3f, n_x=%d, n_y=%d, local=%.3f, running total = %.5f\n",
							t, continuousNewStates[t][0], eps_x, n_x, n_y, locals[t], testSum);
				} else {
					System.out.printf("t=%d: eps_x=%.3f, n_x=%d, n_y=%d, local=%.3f, running total = %.5f\n",
							t, eps_x, n_x, n_y, locals[t], testSum);
				}
			}
		}
		if (debug) {
			avNx /= (double)N;
			avNy /= (double)N;
			System.out.printf("Average n_x=%.3f, Average n_y=%.3f\n", avNx, avNy);
		}
		
		return locals;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public double getLastAverage() {
		return mi;
	}
	
	public int getNumObservations() {
		return continuousData.length;
	}
}
