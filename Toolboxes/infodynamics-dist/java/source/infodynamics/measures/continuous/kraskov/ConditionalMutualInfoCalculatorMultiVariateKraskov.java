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

package infodynamics.measures.continuous.kraskov;

import java.util.Calendar;
import java.util.Random;

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateCommon;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.KdTree;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.NearestNeighbourSearcher;
import infodynamics.utils.UnivariateNearestNeighbourSearcher;

/**
 * <p>Computes the differential conditional mutual information of two multivariate
 *  <code>double[][]</code> sets of observations, conditioned on another
 *  (implementing {@link ConditionalMutualInfoCalculatorMultiVariate}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see references below).
 *  The implementation is made using fast-neighbour searches with an
 *  underlying k-d tree algorithm.
 *  This is an abstract class, building on the common code base in
 *  {@link ConditionalMutualInfoMultiVariateCommon},
 *  to gather common functionality between the two
 *  algorithms defined by Kraskov et al.
 *  Two child classes {@link ConditionalMutualInfoCalculatorMultiVariateKraskov1} and
 *  {@link ConditionalMutualInfoCalculatorMultiVariateKraskov2} then
 *  actually implement the two KSG algorithms.</p>
 * 
 * <p>Crucially, the calculation is performed by examining
 * neighbours in the full joint space (as specified by Frenzel and Pompe)
 * rather than two MI calculators.</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link ConditionalMutualInfoCalculatorMultiVariate},
 * with:
 * <ul>
 * 	<li>For constructors see the child classes.</li>
 *  <li>Further properties are defined in {@link #setProperty(String, String)}.</li>
 *  <li>Computed values are in <b>nats</b>, not bits!</li>
 *  </ul>
 * </p>
 * 
 * <p>Finally, note that {@link Cloneable} is implemented allowing clone()
 *  to produce only an automatic shallow copy, which is fine
 *  for the statistical significance calculation it is intended for
 *  (none of the array 
 *  data will be changed there).</p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Frenzel and Pompe, <a href="http://dx.doi.org/10.1103/physrevlett.99.204101">
 * 	"Partial Mutual Information for Coupling Analysis of Multivariate Time Series"</a>,
 * 	Physical Review Letters, <b>99</b>, p. 204101+ (2007).</li>
 * 	<li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @author Ipek Özdemir
 */
public abstract class ConditionalMutualInfoCalculatorMultiVariateKraskov 
	extends ConditionalMutualInfoMultiVariateCommon
	implements Cloneable { // See comments on clonability above
	
	/**
	 * we compute distances to the kth neighbour in the joint space
	 */
	protected int k = 4;
		
	/**
	 * The norm type in use (see {@link #PROP_NORM_TYPE})
	 */
	protected int normType = EuclideanUtils.NORM_MAX_NORM;
		
	/**
	 * Property name for the number of K nearest neighbours used in
	 * the KSG algorithm in the full joint space (default 4).
	 */
	public final static String PROP_K = "k";
	/**
	 * Property name for what type of norm to use between data points
	 *  for each marginal variable -- Options are defined by 
	 *  {@link KdTree#setNormType(String)} and the
	 *  default is {@link EuclideanUtils#NORM_MAX_NORM}.
	 */
	public final static String PROP_NORM_TYPE = "NORM_TYPE";
	/**
	 * Property name for an amount of random Gaussian noise to be
	 *  added to the data (default is 0).
	 */
	public static final String PROP_ADD_NOISE = "NOISE_LEVEL_TO_ADD";
	/**
	 * Property name for a dynamics exclusion time window 
	 * otherwise known as Theiler window (see Kantz and Schreiber).
	 * Default is 0 which means no dynamic exclusion window.
	 */
	public static final String PROP_DYN_CORR_EXCL_TIME = "DYN_CORR_EXCL";	
	/**
	 * Property name for the number of parallel threads to use in the
	 *  computation (default is to use all available)
	 */
	public static final String PROP_NUM_THREADS = "NUM_THREADS";
	/**
	 * Valid property value for {@link #PROP_NUM_THREADS} to indicate
	 *  that all available processors should be used. 
	 */
	public static final String USE_ALL_THREADS = "USE_ALL";

	/**
	 * Whether to add an amount of random noise to the incoming data
	 */
	protected boolean addNoise = false;
	/**
	 * Amount of random Gaussian noise to add to the incoming data
	 */
	protected double noiseLevel = 0.0;
	/**
	 * Whether we use dynamic correlation exclusion
	 */
	protected boolean dynCorrExcl = false;
	/**
	 * Size of dynamic correlation exclusion window.
	 */
	protected int dynCorrExclTime = 0;
	/**
	 * Number of parallel threads to use in the computation;
	 *  defaults to use all available.
	 */
	protected int numThreads = Runtime.getRuntime().availableProcessors();
	/**
	 * Private variable to record which KSG algorithm number
	 *  this instance is implementing
	 */
	protected boolean isAlgorithm1 = false;

	/**
	 * protected k-d tree data structure (for fast nearest neighbour searches)
	 *  representing the joint source-dest space
	 */
	protected KdTree kdTreeJoint;
	/**
	 * protected k-d tree data structure (for fast nearest neighbour searches)
	 *  representing the (var1,conditional) space
	 */
	protected KdTree kdTreeVar1Conditional;
	/**
	 * protected univariate neighbour searcher data structure (for fast nearest neighbour searches)
	 *  representing the (var1) space; used only if var1 is univariate
	 */
	protected UnivariateNearestNeighbourSearcher uniNNSearcherVar1;
	/**
	 * protected k-d tree data structure (for fast nearest neighbour searches)
	 *  representing the (var2,conditional) space
	 */
	protected KdTree kdTreeVar2Conditional;
	/**
	 * protected univariate neighbour searcher data structure (for fast nearest neighbour searches)
	 *  representing the (var2) space; used only if var2 is univariate
	 */
	protected UnivariateNearestNeighbourSearcher uniNNSearcherVar2;
	/**
	 * protected data structure (for fast nearest neighbour searches)
	 *  representing the conditional space.
	 *  Could be implemented by either a k-d tree or sorted array, depending
	 *  on whether the conditional is multi-variate or not.
	 */
	protected NearestNeighbourSearcher nnSearcherConditional;
	
	/**
	 * Constant for digamma(k), with k the number of nearest neighbours selected
	 */
	protected double digammaK;
	
	/**
	 * Construct an instance of the KSG conditional MI calculator
	 */
	public ConditionalMutualInfoCalculatorMultiVariateKraskov() {
		super();
	}

	@Override
	public void initialise(int dimensions1, int dimensions2, int dimensionsCond) {
		kdTreeJoint = null;
		kdTreeVar1Conditional = null;
		kdTreeVar2Conditional = null;
		nnSearcherConditional = null;
		uniNNSearcherVar1 = null;
		uniNNSearcherVar2 = null;
		super.initialise(dimensions1, dimensions2, dimensionsCond);
	}

	/**
	 * Sets properties for the KSG conditional MI calculator.
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 *  
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 *  <li>{@link #PROP_K} -- number of k nearest neighbours to use in joint kernel space
	 *      in the KSG algorithm (default is 4).</li>
	 *  <li>{@link #PROP_NORM_TYPE}</li> -- normalization type to apply to 
	 *      working out the norms between the points in each marginal space.
	 *      Options are defined by {@link KdTree#setNormType(String)} -
	 *      default is {@link EuclideanUtils#NORM_MAX_NORM}.
	 *  <li>{@link #DYN_CORR_EXCL_TIME_NAME} -- a dynamics exclusion time window,
	 *      also known as Theiler window (see Kantz and Schreiber);
	 *      default is 0 which means no dynamic exclusion window.</li>
	 *  <li>{@link #PROP_ADD_NOISE} -- a standard deviation for an amount of
	 *  	random Gaussian noise to add to
	 *      each variable, to avoid having neighbourhoods with artificially
	 *      large counts. The amount is added in after any normalisation,
	 *      so can be considered as a number of standard deviations of the data.
	 *      (Recommended by Kraskov. MILCA uses 1e-8; but adds in
	 *      a random amount of noise in [0,noiseLevel) ). Default 0.</li>
	 *  <li>{@link #PROP_NUM_THREADS} -- the integer number of parallel threads
	 *  	to use in the computation. Can be passed as a string "USE_ALL"
	 *      to use all available processors on the machine.
	 *      Default is "USE_ALL".
	 *  <li>any valid properties for {@link ConditionalMutualInfoMultiVariateCommon#setProperty(String, String)},
	 *     notably including {@link ConditionalMutualInfoMultiVariateCommon#PROP_NORMALISE}.</li>
	 * </ul>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue) {
		if (propertyName.equalsIgnoreCase(PROP_K)) {
			k = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NORM_TYPE)) {
			normType = KdTree.validateNormType(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_DYN_CORR_EXCL_TIME)) {
			dynCorrExclTime = Integer.parseInt(propertyValue);
			dynCorrExcl = (dynCorrExclTime > 0);
		} else if (propertyName.equalsIgnoreCase(PROP_ADD_NOISE)) {
			addNoise = true;
			noiseLevel = Double.parseDouble(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_NUM_THREADS)) {
			if (propertyValue.equalsIgnoreCase(USE_ALL_THREADS)) {
				numThreads = Runtime.getRuntime().availableProcessors();
			} else { // otherwise the user has passed in an integer:
				numThreads = Integer.parseInt(propertyValue);
			}
		} else {
			// Assume this is a property for the common parent class
			super.setProperty(propertyName, propertyValue);
		}
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateCommon#finaliseAddObservations()
	 */
	@Override
	public void finaliseAddObservations() throws Exception {
		// Allow the parent to generate the data for us first
		super.finaliseAddObservations();
		
		if (dynCorrExcl && addedMoreThanOneObservationSet) {
			// We have not properly implemented dynamic correlation exclusion for
			//  multiple observation sets, so throw an error
			throw new RuntimeException("Addition of multiple observation sets is not currently " +
					"supported with property " + PROP_DYN_CORR_EXCL_TIME + " set");
		}
		
		if (totalObservations <= k + 2*dynCorrExclTime) {
			throw new Exception("There are less observations provided (" +
					totalObservations +
					") than required for the number of nearest neighbours parameter (" +
					k + ") and any dynamic correlation exclusion (" + dynCorrExclTime + ")");
		}
		
		if (addNoise) {
			Random random = new Random();
			// Add Gaussian noise of std dev noiseLevel to the data
			for (int r = 0; r < var1Observations.length; r++) {
				for (int c = 0; c < dimensionsVar1; c++) {
					var1Observations[r][c] +=
							random.nextGaussian()*noiseLevel;
				}
				for (int c = 0; c < dimensionsVar2; c++) {
					var2Observations[r][c] +=
							random.nextGaussian()*noiseLevel;
				}
				for (int c = 0; c < dimensionsCond; c++) {
					condObservations[r][c] +=
							random.nextGaussian()*noiseLevel;
				}
			}
		}

		// Set the constants:
		digammaK = MathsUtils.digamma(k);
	}

	/**
	 * {@inheritDoc} 
	 * 
	 * @return the average conditional MI in nats (not bits!)
	 */
	@Override
	public double computeAverageLocalOfObservations() throws Exception {
		// Compute the conditional MI
		double startTime = Calendar.getInstance().getTimeInMillis();
		lastAverage = computeFromObservations(false)[0];
		condMiComputed = true;
		if (debug) {
			Calendar rightNow2 = Calendar.getInstance();
			long endTime = rightNow2.getTimeInMillis();
			System.out.println("Calculation time: " + ((endTime - startTime)/1000.0) + " sec" );
		}
		return lastAverage;
	}

	/**
	 * {@inheritDoc} 
	 * 
	 * If {@code reordering} is null, it is assumed there is no reordering of
	 *  the given variable.
	 *  
	 * @return the conditional MI under the new ordering, in nats (not bits!).
	 */
	@Override
	public double computeAverageLocalOfObservations(int variableToReorder, 
			int[] reordering) throws Exception {
		
		if (reordering == null) {
			return computeAverageLocalOfObservations();
		}
		double[][] originalData;
		KdTree originalKdTreeJoint = kdTreeJoint;
		kdTreeJoint = null; // So that it is rebuilt for the new ordering
		KdTree originalKdTreeVar1Conditional = kdTreeVar1Conditional;
		UnivariateNearestNeighbourSearcher originalUniNNSearcherVar1 = uniNNSearcherVar1;
		KdTree originalKdTreeVar2Conditional = kdTreeVar2Conditional;
		UnivariateNearestNeighbourSearcher originalUniNNSearcherVar2 = uniNNSearcherVar2;
		if (variableToReorder == 1) {
			originalData = var1Observations;
			kdTreeVar1Conditional = null; // So that it is rebuilt for the new ordering if required
			uniNNSearcherVar1 = null; // So that it is rebuilt for the new ordering if required
			var1Observations = MatrixUtils.extractSelectedTimePointsReusingArrays(originalData, reordering);
		} else {
			originalData = var2Observations;
			kdTreeVar2Conditional = null; // So that it is rebuilt for the new ordering
			uniNNSearcherVar2 = null; // So that it is rebuilt for the new ordering if required
			var2Observations = MatrixUtils.extractSelectedTimePointsReusingArrays(originalData, reordering);
		}
		// Compute the conditional MI
		double newCondMI = computeFromObservations(false)[0];
		// restore original data
		kdTreeJoint = originalKdTreeJoint;
		if (variableToReorder == 1) {
			var1Observations = originalData;
			kdTreeVar1Conditional = originalKdTreeVar1Conditional;
			uniNNSearcherVar1 = originalUniNNSearcherVar1;
		} else {
			var2Observations = originalData;
			kdTreeVar2Conditional = originalKdTreeVar2Conditional;
			uniNNSearcherVar2 = originalUniNNSearcherVar2;
		}
		return newCondMI;
	}

	@Override
	public double[] computeLocalOfPreviousObservations() throws Exception {
		double[] localValues = computeFromObservations(true);
		lastAverage = MatrixUtils.mean(localValues);
		condMiComputed = true;
		return localValues;
	}

	/**
	 * @returns the series of local conditional MI values in nats, not bits.
	 */
	@Override
	public double[] computeLocalUsingPreviousObservations(double[][] states1,
			double[][] states2, double[][] condStates) throws Exception {
		// If implemented, will need to incorporate any normalisation here
		//  (normalising the incoming data the same way the previously
		//   supplied observations were normalised).
		throw new Exception("Local method not implemented yet");
	}
	
	/**
	 * This protected method handles the multiple threads which
	 *  computes either the average or local conditional MI (over parts of the total
	 *  observations), computing the x, y and z 
	 *  distances between all tuples in time.
	 * 
	 * <p>The method returns:<ol>
	 *  <li>for (returnLocals == false), an array of size 1,
	 *      containing the average conditional MI </li>
	 *  <li>for local conditional MIs (returnLocals == true), the array of local conditional MI values</li>
	 *  </ol>
	 * 
	 * @param returnLocals whether to return an array or local values, or else
	 *  sums of these values
	 * @return either the average conditional MI, or array of local conditional MI value, in nats not bits
	 * @throws Exception
	 */
	protected double[] computeFromObservations(boolean returnLocals) throws Exception {
		int N = var1Observations.length; // number of observations
		
		double[] returnValues = null;
		
		// We need to construct the k-d trees for use by the child
		//  classes. We check each tree for existence separately
		//  some can be used across original and surrogate data
		// TODO can parallelise these -- best done within the kdTree --
		//  though it's unclear if there's much point given that
		//  the tree construction itself afterwards can't really be well parallelised.
		if (kdTreeJoint == null) {
			kdTreeJoint = new KdTree(
					new int[] {dimensionsVar1, dimensionsVar2, dimensionsCond},
					new double[][][] {var1Observations, var2Observations, condObservations});
			kdTreeJoint.setNormType(normType);
		}
		if (dimensionsVar1 > 1) {
			if (kdTreeVar1Conditional == null) {
				kdTreeVar1Conditional = new KdTree(
						new int[] {dimensionsVar1, dimensionsCond},
						new double[][][] {var1Observations, condObservations});
				kdTreeVar1Conditional.setNormType(normType);
			} 
		} else { // Univariate variable 1, so we'll search its space alone as this is faster
			if (uniNNSearcherVar1 == null) {
				uniNNSearcherVar1 = new UnivariateNearestNeighbourSearcher(var1Observations);
			}
		}
		if (dimensionsVar2 > 1) {
			if (kdTreeVar2Conditional == null) {
				kdTreeVar2Conditional = new KdTree(
						new int[] {dimensionsVar2, dimensionsCond},
						new double[][][] {var2Observations, condObservations});
				kdTreeVar2Conditional.setNormType(normType);
			}
		} else { // Univariate variable 2, so we'll search its space alone as this is faster
			if (uniNNSearcherVar2 == null) {
				uniNNSearcherVar2 = new UnivariateNearestNeighbourSearcher(var2Observations);
			}
		}
		if (nnSearcherConditional == null) {
			nnSearcherConditional = NearestNeighbourSearcher.create(condObservations);
			nnSearcherConditional.setNormType(normType);
		}

		if (numThreads == 1) {
			// Single-threaded implementation:
			returnValues = partialComputeFromObservations(0, N, returnLocals);
			
		} else {
			// We're going multithreaded:
			if (returnLocals) {
				// We're computing local conditional MI
				returnValues = new double[N];
			} else {
				// We're computing average conditional MI
				returnValues = new double[CondMiKraskovThreadRunner.RETURN_ARRAY_LENGTH];
			}
			
			// Distribute the observations to the threads for the parallel processing
			int lTimesteps = N / numThreads; // each thread gets the same amount of data
			int res = N % numThreads; // the first thread gets the residual data
			if (debug) {
				System.out.printf("Computing Kraskov conditional MI with %d threads (%d timesteps each, plus %d residual)\n",
						numThreads, lTimesteps, res);
			}
			Thread[] tCalculators = new Thread[numThreads];
			CondMiKraskovThreadRunner[] runners = new CondMiKraskovThreadRunner[numThreads];
			for (int t = 0; t < numThreads; t++) {
				int startTime = (t == 0) ? 0 : lTimesteps * t + res;
				int numTimesteps = (t == 0) ? lTimesteps + res : lTimesteps;
				if (debug) {
					System.out.println(t + ".Thread: from " + startTime +
							" to " + (startTime + numTimesteps)); // Trace Message
				}
				runners[t] = new CondMiKraskovThreadRunner(this, startTime, numTimesteps, returnLocals);
				tCalculators[t] = new Thread(runners[t]);
				tCalculators[t].start();
			}
			
			// Here, we should wait for the termination of the all threads
			//  and collect their results
			for (int t = 0; t < numThreads; t++) {
				if (tCalculators[t] != null) { // TODO Ipek: can you comment on why we're checking for null here?
					tCalculators[t].join(); 
				}
				// Now we add in the data from this completed thread:
				if (returnLocals) {
					// We're computing local MI; copy these local values
					//  into the full array of locals
					System.arraycopy(runners[t].getReturnValues(), 0, 
							returnValues, runners[t].myStartTimePoint, runners[t].numberOfTimePoints);
				} else {
					// We're computing the average MI, keep the running sums of digammas and counts
					MatrixUtils.addInPlace(returnValues, runners[t].getReturnValues());
				}
			}
		}
		
		// Finalise the results:
		if (returnLocals) {
			return returnValues;
		} else {
			// Average out the components for the final equation(s) and for debugging:
			double averageDiGammas = returnValues[CondMiKraskovThreadRunner.INDEX_SUM_DIGAMMAS] / (double) N;
			double avNxz = returnValues[CondMiKraskovThreadRunner.INDEX_SUM_NXZ] / (double) N;
			double avNyz = returnValues[CondMiKraskovThreadRunner.INDEX_SUM_NYZ] / (double) N;
			double avNz = returnValues[CondMiKraskovThreadRunner.INDEX_SUM_NZ] / (double) N;
			if (debug) {
				System.out.printf("<n_xz>=%.3f, <n_yz>=%.3f, <n_z>=%.3f\n",
						avNxz, avNyz, avNz);
			}
			if (this.isAlgorithm1) {
				// Algorithm 1:
				if (debug) {
					System.out.printf("Av = digamma(k)=%.3f + <digammas>=%.3f = %.3f \n",
							MathsUtils.digamma(k), averageDiGammas, MathsUtils.digamma(k) + averageDiGammas);
				}
				double[] result = new double[1];
				result[0] = MathsUtils.digamma(k) + averageDiGammas;
				return result;
			} else {
				// Algorithm 2:
				// We also retrieve the sums of inverses for debugging purposes:
				double averageInverseCountInJointXZ =
						returnValues[CondMiKraskovThreadRunner.INDEX_SUM_INV_NXZ] / (double) N;
				double averageInverseCountInJointYZ =
						returnValues[CondMiKraskovThreadRunner.INDEX_SUM_INV_NYZ] / (double) N;
				double averageMeasure = MathsUtils.digamma(k) - (2.0 / (double)k) + averageDiGammas +
						averageInverseCountInJointXZ + averageInverseCountInJointYZ;
				if (debug) {
					System.out.printf("Av = digamma(k)=%.3f + <digammas>=%.3f +<inverses>=%.3f - 2/k=%.3f  = %.3f" +
							" (<1/n_yz>=%.3f, <1/n_xz>=%.3f)\n",
							MathsUtils.digamma(k), averageDiGammas,
							averageInverseCountInJointXZ + averageInverseCountInJointYZ,
							(2.0 / (double)k), averageMeasure,
							averageInverseCountInJointYZ, averageInverseCountInJointXZ);
				}
				double[] result = new double[1];
				result[0] = averageMeasure;
				return result;
			}			
		}
	}
	
	/**
	 * Protected method to be used internally for threaded implementations.
	 * This method implements the guts of each Kraskov algorithm, computing the number of 
	 *  nearest neighbours in each dimension for a sub-set of the data points.
	 *  It is intended to be called by one thread to work on that specific
	 *  sub-set of the data.
	 * 
	 * <p>The method returns:<ol>
	 *  <li>for average conditional MIs (returnLocals == false), the relevant sums of digamma(n_{xz}), digamma(n_{yz})
	 *     and digamma(n_z)
	 *     for a partial set of the observations</li>
	 *  <li>for local conditional MIs (returnLocals == true), the array of local conditional MI values</li>
	 *  </ol>
	 * 
	 * @param startTimePoint start time for the partial set we examine
	 * @param numTimePoints number of time points (including startTimePoint to examine)
	 * @param returnLocals whether to return an array or local values, or else
	 *  sums of these values
	 * @return an array of the relevant sum of digamma(n_xz+1) and digamma(n_yz+1) and digamma(n_z), then
	 *  sum of n_xz, n_yz, n_z and for algorithm 2 only, sum of 1/n_xz and 1/n_yz
	 *  (these latter five are only for debugging purposes).
	 * @throws Exception
	 */
	protected abstract double[] partialComputeFromObservations(
			int startTimePoint, int numTimePoints, boolean returnLocals) throws Exception;

	/**
	 * Private class to handle multi-threading of the Kraskov algorithms.
	 * Each instance calls partialComputeFromObservations()
	 * to compute nearest neighbours for a part of the data.
	 * 
	 * 
	 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
	 * <a href="http://lizier.me/joseph/">www</a>)
	 * @author Ipek Özdemir
	 */
	private class CondMiKraskovThreadRunner implements Runnable {
		protected ConditionalMutualInfoCalculatorMultiVariateKraskov condMiCalc;
		protected int myStartTimePoint;
		protected int numberOfTimePoints;
		protected boolean computeLocals;
		
		protected double[] returnValues = null;
		protected Exception problem = null;
		
		public static final int INDEX_SUM_DIGAMMAS = 0;
		public static final int INDEX_SUM_NXZ = 1;
		public static final int INDEX_SUM_NYZ = 2;
		public static final int INDEX_SUM_NZ = 3;
		public static final int INDEX_SUM_INV_NXZ = 4; // Only used for algorithm 2
		public static final int INDEX_SUM_INV_NYZ = 5; // Only used for algorithm 2
		public static final int RETURN_ARRAY_LENGTH = 6;
		
		public CondMiKraskovThreadRunner(
				ConditionalMutualInfoCalculatorMultiVariateKraskov condMiCalc,
				int myStartTimePoint, int numberOfTimePoints,
				boolean computeLocals) {
			this.condMiCalc = condMiCalc;
			this.myStartTimePoint = myStartTimePoint;
			this.numberOfTimePoints = numberOfTimePoints;
			this.computeLocals = computeLocals;
		}
		
		/**
		 * Return the values from this part of the data,
		 *  or throw any exception that was encountered by the 
		 *  thread.
		 * 
		 * @return an exception previously encountered by this thread.
		 * @throws Exception
		 */
		public double[] getReturnValues() throws Exception {
			if (problem != null) {
				throw problem;
			}
			return returnValues;
		}
		
		/**
		 * Start the thread for the given parameters
		 */
		public void run() {
			try {
				returnValues = condMiCalc.partialComputeFromObservations(myStartTimePoint, numberOfTimePoints, computeLocals);
			} catch (Exception e) {
				// Store the exception for later retrieval
				problem = e;
				return;
			}
		}
	}
	// end class MiKraskovThreadRunner

	// Note: no extra implementation of clone provided; we're simply
	//  allowing clone() to produce a shallow copy, which is fine
	//  for the statistical significance calculation (none of the array
	//  data will be changed there.
	//
	// public ConditionalMutualInfoCalculatorMultiVariateKraskov clone() {
	//	return this;
	// }
}
