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

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.measures.mixed.ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSourceCommon;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Compute the Conditional Mutual Information between a discrete variable and a 
 *  vector of continuous variables, conditioned on another vector of continuous variables
 *  using the Kraskov estimation method.</p>
 * <p>Uses Kraskov method type 2, since type 1 only looks at points with
 * distances strictly less than the kth variable, which won't work for one marginal
 * being discrete.</p>
 * <p>The actual equation which should be used for type 2 with conditional MI 
 * follows {#link ConditionalMutualInfoCalculatorMultiVariateKraskov2}.
 * </p>
 * 
 * <p>These calculators are <b>EXPERIMENTAL</b> -- not properly tested,
 * and not well documented. The intended calling pattern is similar to
 * {@link ConditionalMutualInfoCalculatorMultiVariate}
 * </p>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class ConditionalMutualInfoCalculatorMultiVariateWithDiscreteKraskov
	extends ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSourceCommon
	implements Cloneable { // See comments on clonability below

	// Multiplier used in hueristic for determining whether to use a linear search
	//  for min kth element or a binary search.
	protected static final double CUTOFF_MULTIPLIER = 1.5;
	
	/**
	 * we compute distances to the kth neighbour
	 */
	protected int k = 4;
	
	protected EuclideanUtils normCalculator;
	// Storage for the norms from each observation to each other one
	protected double[][] xNorms;
	protected double[][] zNorms;
	protected double[][] xzNorms;
	// Keep the norms each time (making reordering very quick)
	//  (Should only be set to false for testing)
	public static boolean tryKeepAllPairsNorms = true;
	public static int MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM = 2000;
	
	public final static String PROP_K = "k";
	public final static String PROP_NORM_TYPE = "NORM_TYPE";	

	public ConditionalMutualInfoCalculatorMultiVariateWithDiscreteKraskov() {
		super();
		normCalculator = new EuclideanUtils(EuclideanUtils.NORM_MAX_NORM);
	}

	/**
	 * Initialise the calculator.
	 * 
	 * @param dimensions number of joint continuous variables
	 * @param base number of discrete states
	 * @param dimensionsCond the number of joint continuous variables
	 *     to condition on
	 */
	public void initialise(int dimensions, int base, int dimensionsCond) {
		super.initialise(dimensions, base, dimensionsCond);

		xNorms = null;
		zNorms = null;
		xzNorms = null;
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
	 *  <li>Any other properties settable in the parent class'
	 *  {@link ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSourceCommon#setProperty(String, String)}</li>
	 * </ul>
	 * 
	 * @param propertyName
	 * @param propertyValue
	 */
	public void setProperty(String propertyName, String propertyValue) {
		if (propertyName.equalsIgnoreCase(PROP_K)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
			normCalculator.setNormToUse(propertyValue);
		} else {
			super.setProperty(propertyName, propertyValue);
		}
	}

	public void finaliseAddObservations() throws Exception {
		super.finaliseAddObservations();
		// Now check that we have at least k observations in each discrete bin, 
		//  or else our Kraskov extension won't make sense:
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
		int N = continuousDataX.length; // number of observations
		
		xNorms = new double[N][N];
		zNorms = new double[N][N];
		xzNorms = new double[N][N];
		for (int t = 0; t < N; t++) {
			// Compute the norms from t to all other time points
			for (int t2 = 0; t2 < N; t2++) {
				if (t2 == t) {
					xNorms[t][t2] = Double.POSITIVE_INFINITY;
					zNorms[t][t2] = Double.POSITIVE_INFINITY;
					xzNorms[t][t2] = Double.POSITIVE_INFINITY;
					continue;
				}
				// Compute norm in the continuous space
				xNorms[t][t2] = normCalculator.norm(continuousDataX[t], continuousDataX[t2]);
				zNorms[t][t2] = normCalculator.norm(conditionedDataZ[t], conditionedDataZ[t2]);
				xzNorms[t][t2] = Math.max(xNorms[t][t2], zNorms[t][t2]);
			}
		}
	}
	
	/**
	 * Compute what the average conditional MI would look like were the second time series reordered
	 *  as per the array of time indices in reordering.
	 * The user should ensure that all values 0..N-1 are represented exactly once in the
	 *  array reordering and that no other values are included here. 
	 * 
	 * @param reordering
	 * @return
	 * @throws Exception
	 */
	public double computeAverageLocalOfObservations(int[] reordering) throws Exception {
		int N = continuousDataX.length; // number of observations
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
		double avNxz = 0;
		double avNyz = 0;
		double avNz = 0;
		double averageInverseCountInJointXZ = 0;
		double averageInverseCountInJointYZ = 0;

		for (int t = 0; t < N; t++) {
			// Compute eps_x and eps_z for this time step:
			//  using max of x and z norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).

			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2][0] = Math.max(xNorms[t][t2], zNorms[t][t2]);
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][1] = t2;
			}
			// Then find the k closest neighbours:
			double eps_x = 0.0;
			double eps_z = 0.0;
			int[] timeStepsOfKthMins = null;
			// just do a linear search for the minimum epsilon value
			timeStepsOfKthMins = MatrixUtils.kMinIndicesSubjectTo(
					jointNorm, 0, k, reorderedDiscreteData, reorderedDiscreteData[t]);
			// and now we have the closest k points.
			// Find eps_{x,y,z} as the maximum x and y and z norms amongst this set:
			for (int j = 0; j < k; j++) {
				int timeStepOfJthPoint = timeStepsOfKthMins[j]; 
				if (xNorms[t][timeStepOfJthPoint] > eps_x) {
					eps_x = xNorms[t][timeStepOfJthPoint];
				}
				if (zNorms[t][timeStepOfJthPoint] > eps_z) {
					eps_z = zNorms[t][timeStepOfJthPoint];
				}
			}

			// Count the number of points whose distances are less
			//  than or equal to eps in each required joint space
			int n_xz = 0;
			int n_yz = 0;
			int n_z = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (zNorms[t][t2] <= eps_z) {
					n_z++;
					if (xNorms[t][t2] <= eps_x) {
						n_xz++;
					}
					if (reorderedDiscreteData[t] == reorderedDiscreteData[t2]) {
						n_yz++;
					}
				}
			}
			avNxz += n_xz;
			avNyz += n_yz;
			avNz += n_z;
			// And take the digamma before adding into the 
			//  average:
			averageDiGammas += MathsUtils.digamma(n_z)
					- MathsUtils.digamma(n_xz) - MathsUtils.digamma(n_yz);
			double invN_xz = 1.0/(double) n_xz;
			averageInverseCountInJointXZ += invN_xz;
			double invN_yz = 1.0 / (double) n_yz;
			averageInverseCountInJointYZ += invN_yz;
		}
		averageDiGammas /= (double) N;
		averageInverseCountInJointXZ /= (double) N;
		averageInverseCountInJointYZ /= (double) N;
		condMi = MathsUtils.digamma(k) - 2.0 / (double) k + 
				averageDiGammas + averageInverseCountInJointXZ +
				averageInverseCountInJointYZ;
		miComputed = true;
		
		if (debug) {
			avNxz /= (double)N;
			avNyz /= (double)N;
			avNz /= (double)N;
			System.out.println(String.format("Average n_xz=%.3f, Average n_yz=%.3f, Average n_z=%.3f",
					avNxz, avNyz, avNz));
			System.out.printf("Av = digamma(k)=%.3f + <digammas>=%.3f +<inverses>=%.3f - 2/k=%.3f = %.3f (<1/n_yz>=%.3f, <1/n_xz>=%.3f)\n",
					MathsUtils.digamma(k), averageDiGammas, 
					averageInverseCountInJointXZ + averageInverseCountInJointYZ,
					2.0 / (double) k,
					condMi, averageInverseCountInJointYZ, averageInverseCountInJointXZ);
		}
		
		return condMi;
	}

	public double computeAverageLocalOfObservations() throws Exception {
		if (!tryKeepAllPairsNorms || (continuousDataX.length > MAX_DATA_SIZE_FOR_KEEP_ALL_PAIRS_NORM)) {
			return computeAverageLocalOfObservationsWhileComputingDistances();
		}
		
		if (xNorms == null) {
			computeNorms();
		}
		int N = continuousDataX.length; // number of observations

		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double averageInverseCountInJointXZ = 0;
		double averageInverseCountInJointYZ = 0;
		double avNxz = 0;
		double avNyz = 0;
		double avNz = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps_x and eps_z for this time step:
			//  using x,z norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).

			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2][0] = Math.max(xNorms[t][t2], zNorms[t][t2]);
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][1] = t2;
			}
			// Then find the k closest neighbours:
			double eps_x = 0.0;
			double eps_z = 0.0;
			int[] timeStepsOfKthMins = null;
			// just do a linear search for the minimum epsilon value
			timeStepsOfKthMins = MatrixUtils.kMinIndicesSubjectTo(
					jointNorm, 0, k, discreteData, discreteData[t]);
			// and now we have the closest k points.
			// Find eps_{x,y,z} as the maximum x and y and z norms amongst this set:
			for (int j = 0; j < k; j++) {
				int timeStepOfJthPoint = timeStepsOfKthMins[j]; 
				if (xNorms[t][timeStepOfJthPoint] > eps_x) {
					eps_x = xNorms[t][timeStepOfJthPoint];
				}
				if (zNorms[t][timeStepOfJthPoint] > eps_z) {
					eps_z = zNorms[t][timeStepOfJthPoint];
				}
			}

			// Count the number of points whose x,y,z distances are less
			//  than or equal to eps (not including this point)
			int n_xz = 0;
			int n_yz = 0;
			int n_z = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (zNorms[t][t2] <= eps_z) {
					n_z++;
					if (xNorms[t][t2] <= eps_x) {
						n_xz++;
					}
					if (discreteData[t] == discreteData[t2]) {
						n_yz++;
					}
				}
			}
			avNxz += n_xz;
			avNyz += n_yz;
			avNz += n_z;
			// And take the digamma before adding into the 
			//  average:
			averageDiGammas += MathsUtils.digamma(n_z)
					- MathsUtils.digamma(n_xz) - MathsUtils.digamma(n_yz);
			double invN_xz = 1.0/(double) n_xz;
			averageInverseCountInJointXZ += invN_xz;
			double invN_yz = 1.0 / (double) n_yz;
			averageInverseCountInJointYZ += invN_yz;
		}
		averageDiGammas /= (double) N;
		averageInverseCountInJointYZ /= (double) N;
		averageInverseCountInJointXZ /= (double) N;
		condMi = MathsUtils.digamma(k) - 2.0 / (double) k +
				averageDiGammas + averageInverseCountInJointYZ +
				averageInverseCountInJointXZ;
		miComputed = true;
		
		if (debug) {
			avNxz /= (double)N;
			avNyz /= (double)N;
			avNz /= (double) N;
			System.out.printf("Average n_xz=%.3f (-> digam=%.3f %.3f), Average n_yz=%.3f (-> digam=%.3f)",
					avNxz, MathsUtils.digamma((int) avNxz), MathsUtils.digamma((int) avNxz - 1), avNyz, MathsUtils.digamma((int) avNyz));
			System.out.printf(", Average n_z=%.3f (-> digam=%.3f)\n", avNz, MathsUtils.digamma((int) avNz));
			System.out.printf("Independent average num in joint box is %.3f\n", (avNxz * avNyz / (double) N));
			System.out.printf("Av = digamma(k)=%.3f + <digammas>=%.3f + <avInverses>=%.3f - 2/k=%.3f = %.3f (<1/n_yz>=%.3f, <1/n_xz>=%.3f)\n",
					MathsUtils.digamma(k), averageDiGammas,
					averageInverseCountInJointYZ + averageInverseCountInJointXZ,
					2.0 / (double) k,
					condMi, averageInverseCountInJointYZ, averageInverseCountInJointXZ);
		}
		
		return condMi;
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
		int N = continuousDataX.length; // number of observations

		// Count the average number of points within eps_x and eps_y
		double averageDiGammas = 0;
		double averageInverseCountInJointYZ = 0;
		double averageInverseCountInJointXZ = 0;
		double avNxz = 0;
		double avNyz = 0;
		double avNz = 0;
		
		for (int t = 0; t < N; t++) {
			// Compute eps_* for this time step:
			//  First get xz norms to all neighbours
			//  (note that norm of point t to itself will be set to infinity).
			double[][] xzNorms = normCalculator.computeNorms(continuousDataX, conditionedDataZ, t);
			double[][] jointNorm = new double[N][2];
			for (int t2 = 0; t2 < N; t2++) {
				jointNorm[t2][0] = Math.max(xzNorms[t2][0], 
						xzNorms[t2][1]);
				// And store the time step for back reference after the 
				//  array is sorted.
				jointNorm[t2][1] = t2;
			}
			// Then find the k closest neighbours:
			double eps_x = 0.0;
			double eps_z = 0.0;
			int[] timeStepsOfKthMins = null;
			// just do a linear search for the minimum epsilon value
			//  subject to the discrete variable value
			timeStepsOfKthMins = MatrixUtils.kMinIndicesSubjectTo(
					jointNorm, 0, k, discreteData, discreteData[t]);
			// and now we have the closest k points.
			// Find eps_{x,z} as the maximum x and z norms amongst this set:
			for (int j = 0; j < k; j++) {
				int timeStepOfJthPoint = timeStepsOfKthMins[j];
				if (xzNorms[timeStepOfJthPoint][0] > eps_x) {
					eps_x = xzNorms[timeStepOfJthPoint][0];
				}
				if (xzNorms[timeStepOfJthPoint][1] > eps_z) {
					eps_z = xzNorms[timeStepOfJthPoint][1];
				}
			}

			// Count the number of points whose distances is less
			//  than or equal to eps
			int n_xz = 0;
			int n_yz = 0;
			int n_z = 0;
			for (int t2 = 0; t2 < N; t2++) {
				if (xzNorms[t2][1] <= eps_z) {
					n_z++;
					if (xzNorms[t2][0] <= eps_x) {
						n_xz++;
					}
					if (discreteData[t] == discreteData[t2]) {
						n_yz++;
					}
				}
			}
			avNxz += n_xz;
			avNyz += n_yz;
			avNz += n_z;
			// And take the digamma before adding into the 
			//  average:
			averageDiGammas += MathsUtils.digamma(n_z)
					- MathsUtils.digamma(n_xz) - MathsUtils.digamma(n_yz);
			double invN_xz = 1.0/(double) n_xz;
			averageInverseCountInJointXZ += invN_xz;
			double invN_yz = 1.0 / (double) n_yz;
			averageInverseCountInJointYZ += invN_yz;
		}
		averageDiGammas /= (double) N;
		averageInverseCountInJointYZ /= (double) N;
		averageInverseCountInJointXZ /= (double) N;
		condMi = MathsUtils.digamma(k) - 2.0 / (double) k +
				averageDiGammas + averageInverseCountInJointYZ +
				averageInverseCountInJointXZ;
		miComputed = true;
		
		if (debug) {
			avNxz /= (double)N;
			avNyz /= (double)N;
			avNz /= (double)N;
			System.out.printf("Average n_xz=%.3f, Average n_yz=%.3f, Average n_z=%.3f\n",
					avNxz, avNyz, avNz);
			System.out.printf("Av = digamma(k)=%.3f + <digammas>=%.3f + <inverses>=%.3f - 2/k=%.3f = %.3f (<1/n_yz>=%.3f, <1/n_xz>=%.3f)\n",
					MathsUtils.digamma(k), averageDiGammas,
					averageInverseCountInJointYZ + averageInverseCountInJointXZ,
					2.0 / (double) k,
					condMi, averageInverseCountInJointYZ, averageInverseCountInJointXZ);
		}
		
		return condMi;
	}

	public double[] computeLocalOfPreviousObservations() throws Exception {
		throw new Exception("Not implemented yet");
	}

	public double[] computeLocalUsingPreviousObservations(
			double[][] contNewStates, int[] discreteNewStates,
			double[][] conditionedNewStates) throws Exception {
		
		if (normalise) {
			contNewStates = MatrixUtils.normaliseIntoNewArray(
					contNewStates, meansX, stdsX);
			conditionedNewStates = MatrixUtils.normaliseIntoNewArray(
					conditionedNewStates, meansZ, stdsZ);
		}
		
		throw new Exception("Not implemented yet");
	}

	// Note: no extra implementation of clone provided; we're simply
	//  allowing clone() to produce a shallow copy, which is find
	//  for the statistical significance calculation (none of the array
	//  data will be changed there.
	//
	// public ConditionalMutualInfoCalculatorMultiVariateKraskov clone() {
	//	return this;
	// }
	
}
