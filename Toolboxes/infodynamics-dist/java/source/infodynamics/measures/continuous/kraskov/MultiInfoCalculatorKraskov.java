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

import infodynamics.measures.continuous.MultiInfoCalculator;
import infodynamics.measures.continuous.MultiInfoCalculatorCommon;
import infodynamics.utils.EuclideanUtils;
import infodynamics.utils.KdTree;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.UnivariateNearestNeighbourSearcher;

import java.util.Calendar;
import java.util.Random;

/**
 * <p>Computes the differential multi-information of a given multivariate set of
 *  observations (implementing {@link MultiInfoCalculator}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see Kraskov et al., below).
 *  This is an abstract class to gather common functionality between the two
 *  algorithms defined by Kraskov et al.
 *  Two child classes {@link MultiInfoCalculatorKraskov1} and
 *  {@link MultiInfoCalculatorKraskov2} then
 *  actually implement the two algorithms in the Kraskov et al. paper</p>
 *    
 * <p>Usage is as per the paradigm outlined for {@link MultiInfoCalculator},
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
 *  data will be changed there).
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 * </ul>
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @author Ipek Özdemir
 */
public abstract class MultiInfoCalculatorKraskov
	extends MultiInfoCalculatorCommon
	implements Cloneable { // See comments on clonability above

	/**
	 * we compute distances to the kth nearest neighbour
	 */
	protected int k = 4;
		
	/**
	 * The norm type in use (see {@link #PROP_NORM_TYPE})
	 */
	protected int normType = EuclideanUtils.NORM_MAX_NORM;
		
	/**
	 * Property name for the number of K nearest neighbours used in
	 * the KSG algorithm (default 4).
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
	 *  representing the joint space
	 */
	protected KdTree kdTreeJoint;
	/**
	 * protected data structures (for fast nearest neighbour searches)
	 *  representing the marginal spaces
	 */
	protected UnivariateNearestNeighbourSearcher[] rangeSearchersInMarginals;
	/**
	 * Constant for digamma(k), with k the number of nearest neighbours selected
	 */
	protected double digammaK;
	/**
	 * Constant for digamma(N), with N the number of samples.
	 */
	protected double digammaN;

	/**
	 * Construct an instance
	 */
	public MultiInfoCalculatorKraskov() {
		super();
	}

	@Override
	public void initialise(int dimensions) {
		// Now call the super class to handle the common variables:
		super.initialise(dimensions);
		kdTreeJoint = null;
		rangeSearchersInMarginals = null;
	}

	/**
	 * Sets properties for the KSG multi-info calculator.
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 *  
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 *  <li>{@link #PROP_K} -- number of k nearest neighbours to use in joint kernel space
	 *      in the KSG algorithm (default is 4).</li>
	 *  <li>{@link #PROP_NORM_TYPE} -- normalization type to apply to 
	 *      working out the norms between the points in each marginal space.
	 * 	    Options are defined by {@link KdTree#setNormType(String)} -
	 * 	    default is {@link EuclideanUtils#NORM_MAX_NORM}.</li>
	 *  <li>{@link #PROP_DYN_CORR_EXCL_TIME} -- a dynamics exclusion time window,
	 *      also known as Theiler window (see Kantz and Schreiber);
	 *      default is 0 which means no dynamic exclusion window.</li>
	 *  <li>{@link #PROP_ADD_NOISE} -- a standard deviation for an amount of
	 *  	random Gaussian noise to add to
	 *      each variable, to avoid having neighbourhoods with artificially
	 *      large counts. The amount is added in after any normalisation,
	 *      so can be considered as a number of standard deviations of the data.
	 *      (Recommended by Kraskov. MILCA uses 1e-8; but adds in
	 *      a random amount of noise in [0,noiseLevel) ). Default 0.</li>
	 * </ul>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		boolean propertySet = true;
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
			// No property was set here
			propertySet = false;
			// try the superclass:
			super.setProperty(propertyName, propertyValue);
		}
		if (debug && propertySet) {
			System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	/**
	 * Set observations from two separate time series: join the rows at each time step
	 *  together to make a joint vector, then effectively call {@link #setObservations(double[][])}
	 * 
	 * @param observations1 observations for first few variables
	 * @param observations2 observations for the other variables
	 * @see #setObservations(double[][])
	 */
	public void setObservations(double[][] observations1, double[][] observations2) throws Exception {
		if ((observations1 == null) || (observations1[0].length == 0) ||
				(observations2 == null) || (observations2[0].length == 0)) {
			throw new Exception("Computing MI with a null set of data");
		}
		if (observations1.length != observations2.length) {
			throw new Exception("Length of the time series to be joined to not match");
		}
		if (observations1[0].length + observations2[0].length != dimensions) {
			throw new Exception("Incorrect number of dimensions " + 
					(observations1[0].length + observations2[0].length) +
					" in supplied observations (expected " + dimensions + ")");
		}
		totalObservations = observations1.length;
		observations = new double[totalObservations][dimensions];
		for (int t = 0; t < totalObservations; t++) {
			int v = 0;
			for (int i = 0; i < observations1[t].length; i++) {
				observations[t][v++] = observations1[t][i];
			}
			for (int i = 0; i < observations2[t].length; i++) {
				observations[t][v++] = observations2[t][i];
			}
		}
		// Normalise the data if required
		if (normalise) {
			// We can overwrite these since they're already
			//  a copy of the users' data.
			MatrixUtils.normalise(observations);
		}
		
		if (addNoise) {
			Random random = new Random();
			// Add Gaussian noise of std dev noiseLevel to the data
			for (int r = 0; r < observations.length; r++) {
				for (int c = 0; c < dimensions; c++) {
					observations[r][c] += random.nextGaussian()*noiseLevel;
				}
			}
		}
	}
	
	@Override
	public void finaliseAddObservations() throws Exception {
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
		
		// Normalise the data if required -- common class doesn't do this
		//  since the Kernel estimator needs to do it internally
		if (normalise) {
			// normalise into new array so we don't overwrite
			//  users' original data
			observations = MatrixUtils.normaliseIntoNewArray(observations);
		}
		
		if (addNoise) {
			Random random = new Random();
			// Add Gaussian noise of std dev noiseLevel to the data
			for (int r = 0; r < observations.length; r++) {
				for (int c = 0; c < dimensions; c++) {
					observations[r][c] += random.nextGaussian()*noiseLevel;
				}
			}
		}

		// Set the constants:
		digammaK = MathsUtils.digamma(k);
		digammaN = MathsUtils.digamma(totalObservations);
	}

	/**
	 * {@inheritDoc} 
	 * 
	 * @return the average multi-info in nats (not bits!)
	 */
	@Override
	public double computeAverageLocalOfObservations() throws Exception {
		// Compute the MI
		double startTime = Calendar.getInstance().getTimeInMillis();
		lastAverage = computeFromObservations(false)[0];
		miComputed = true;
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
	 * @return the "time-series" of local multi-info values in nats (not bits!)
	 * @throws Exception
	 */
	@Override
	public double[] computeLocalOfPreviousObservations() throws Exception {
		double[] localValues = computeFromObservations(true);
		lastAverage = MatrixUtils.mean(localValues);
		miComputed = true;
		return localValues;
	}
	
	/**
	 * This method, specified in {@link MultiInfoCalculator}
	 * is not implemented yet here.
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] states) throws Exception {
		// TODO If this is implemented, will need to normalise the incoming
		//  observations the same way that previously supplied ones were
		//  normalised (if they were normalised, that is)
		throw new Exception("Local method not implemented yet");
	}

	/**
	 * This protected method handles the multiple threads which
	 *  computes either the average or local multi-info (over parts of the total
	 *  observations), computing the
	 *  distances between all tuples in time.
	 * 
	 * <p>The method returns:<ol>
	 *  <li>for (returnLocals == false), an array of size 1,
	 *      containing the average multi-info </li>
	 *  <li>for local multi-infos (returnLocals == true), the array of local
	 *      multi-info values</li>
	 *  </ol>
	 * 
	 * @param returnLocals whether to return an array or local values, or else
	 *  sums of these values
	 * @return either the average multi-info, or array of local multi-info value,
	 * 	in nats not bits
	 * @throws Exception
	 */
	protected double[] computeFromObservations(boolean returnLocals) throws Exception {
		
		double[] returnValues = null;
		
		// We need to construct the k-d trees for use by the child
		//  classes. We check each tree for existence separately
		//  since source can be used across original and surrogate data
		// TODO can parallelise these -- best done within the kdTree --
		//  though it's unclear if there's much point given that
		//  the tree construction itself afterwards can't really be well parallelised.
		double[][][] separateMarginals = null;
		int[] dimensionsArray = null;
		if (kdTreeJoint == null) {
			// We need to pull out 2D time series (of only one variable)
			//  for each marginal variable here
			separateMarginals = new double[dimensions][][];
			dimensionsArray = new int[dimensions];
			for (int d = 0; d < dimensions; d++) {
				separateMarginals[d] =
						MatrixUtils.selectColumns(observations, new int[] {d});
				dimensionsArray[d] = 1;
			}
			kdTreeJoint = new KdTree(dimensionsArray, separateMarginals);
			kdTreeJoint.setNormType(normType);
		}
		if (rangeSearchersInMarginals == null) {
			rangeSearchersInMarginals = new UnivariateNearestNeighbourSearcher[dimensions];
			for (int d = 0; d < dimensions; d++) {
				rangeSearchersInMarginals[d] =
						new UnivariateNearestNeighbourSearcher(
								MatrixUtils.selectColumn(observations, d));
				rangeSearchersInMarginals[d].setNormType(normType);
			}
		}
		
		if (numThreads == 1) {
			// Single-threaded implementation:
			returnValues = partialComputeFromObservations(0, totalObservations, returnLocals);
			
		} else {
			// We're going multithreaded:
			if (returnLocals) {
				// We're computing local MI
				returnValues = new double[totalObservations];
			} else {
				// We're computing average MI
				returnValues = new double[1 + dimensions];
			}
			
			// Distribute the observations to the threads for the parallel processing
			int lTimesteps = totalObservations / numThreads; // each thread gets the same amount of data
			int res = totalObservations % numThreads; // the first thread gets the residual data
			if (debug) {
				System.out.printf("Computing Kraskov Multi-Info with %d threads (%d timesteps each, plus %d residual)\n",
						numThreads, lTimesteps, res);
			}
			Thread[] tCalculators = new Thread[numThreads];
			MultiInfoKraskovThreadRunner[] runners = new MultiInfoKraskovThreadRunner[numThreads];
			for (int t = 0; t < numThreads; t++) {
				int startTime = (t == 0) ? 0 : lTimesteps * t + res;
				int numTimesteps = (t == 0) ? lTimesteps + res : lTimesteps;
				if (debug) {
					System.out.println(t + ".Thread: from " + startTime +
							" to " + (startTime + numTimesteps)); // Trace Message
				}
				runners[t] = new MultiInfoKraskovThreadRunner(this, startTime, numTimesteps, returnLocals);
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
					// We're computing local multi-info; copy these local values
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
			// Compute the average number of points within eps_x and eps_y
			double averageDiGammas = returnValues[MultiInfoKraskovThreadRunner.INDEX_SUM_DIGAMMAS] / (double) totalObservations;
			double[] avNMarginals = new double[dimensions];
			for (int d = 0; d < dimensions; d++) {
				avNMarginals[d] = returnValues[1 + d]  / (double) totalObservations;
				if (debug) {
					System.out.printf("Average n_%d=%.3f, ", d, avNMarginals[d]);
				}
			}
			if (debug) {
				System.out.println();
			}
			
			// Finalise the average result, depending on which algorithm we are implementing:
			if (isAlgorithm1) {
				return new double[] { digammaK - averageDiGammas + (double) (dimensions - 1) * digammaN };
			} else {
				return new double[] { digammaK - ((double) (dimensions - 1)  / (double)k) - averageDiGammas +
									(double) (dimensions - 1) * digammaN };
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
	 *  <li>for average Multi-infos (returnLocals == false), the relevant sums of
	 *  	digamma(n_x+1) in each marginal
	 *     for a partial set of the observations</li>
	 *  <li>for local MIs (returnLocals == true), the array of local MI values</li>
	 *  </ol>
	 * 
	 * @param startTimePoint start time for the partial set we examine
	 * @param numTimePoints number of time points (including startTimePoint to examine)
	 * @param returnLocals whether to return an array or local values, or else
	 *  sums of these values
	 * @return an array of sum of digamma(n_x+1) for each marginal x, then
	 *  sum of n_x for each marginal x (these latter ones are for debugging purposes).
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
	private class MultiInfoKraskovThreadRunner implements Runnable {
		protected MultiInfoCalculatorKraskov miCalc;
		protected int myStartTimePoint;
		protected int numberOfTimePoints;
		protected boolean computeLocals;
		
		protected double[] returnValues = null;
		protected Exception problem = null;
		
		public static final int INDEX_SUM_DIGAMMAS = 0;

		public MultiInfoKraskovThreadRunner(
				MultiInfoCalculatorKraskov miCalc,
				int myStartTimePoint, int numberOfTimePoints,
				boolean computeLocals) {
			this.miCalc = miCalc;
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
				returnValues = miCalc.partialComputeFromObservations(
						myStartTimePoint, numberOfTimePoints, computeLocals);
			} catch (Exception e) {
				// Store the exception for later retrieval
				problem = e;
				return;
			}
		}
	}
	// end class MultiInfoKraskovThreadRunner
}
