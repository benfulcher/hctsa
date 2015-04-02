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

package infodynamics.measures.continuous.kernel;

import infodynamics.utils.MatrixUtils;

import java.util.Arrays;
import java.util.Vector;
import java.util.Hashtable;


/**
 * <p>Class to maintain probability distribution function for
 *  a multivariate set, using kernel estimates.</p>
 * 
 * <p>
 *  For more details on kernel estimation for computing probability distribution functions,
 *  see Kantz and Schreiber (below).
 * </p>
 * 
 * <p>
 * TODO More thoroughly check the Javadocs here
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>H. Kantz and T. Schreiber, "Nonlinear Time Series Analysis"
 *  (Cambridge University Press, Cambridge, MA, 1997).</li>
 * </ul>
 *
 * @see KernelEstimatorUniVariate
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class KernelEstimatorMultiVariate implements Cloneable {

	protected double[] suppliedKernelWidths = null;
	protected double[] kernelWidthsInUse = null;
	protected int dimensions = 1;
	private int[] bins = null;
	private boolean usingIntegerIndexBins = true;
	private boolean useBins = true;
	private double[] mins = null;
	private int[] multipliers = null;
	private int totalObservations = 0;
	private Vector<TimeStampedObservation>[] observations = null;
	private Hashtable<IntArray, Vector<TimeStampedObservation>> observationsByHash = null;
	protected double[][] rawData = null;
	private boolean debug;
	
	// Overriding classes should set this to true if they want a callback for
	//  every point that is found to be correlated in time;
	protected boolean makeCorrelatedPointAddedCallback = false;
	
	protected boolean normalise = true;
	
	private boolean excludeDynamicCorrelations = false;
	private int timeProximityForDynamicCorrelationExclusion = 100;
	
	// Force the kernel estimator to compare each point to every other rather
	//  than use fast bin techniques
	private boolean forceCompareToAll = false;
	
	private final static int MAX_NUMBER_INTEGER_BINS = 1000000;
	
	/**
	 * 
	 * Private class to store a time-stamped data point.
	 * This allows us to eliminate dynamic correlations later.
	 * 
	 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
	 * <a href="http://lizier.me/joseph/">www</a>)
	 */
	private class TimeStampedObservation {
		public int timeStep;
		public double[] observation;
		
		TimeStampedObservation(int time, double[] data) {
			timeStep = time;
			observation = data;
		}
	}
	
	/**
	 * Wrapper class for an integer array so we can hash it properly
	 * 
	 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
	 * <a href="http://lizier.me/joseph/">www</a>)
	 */
	private class IntArray {
		public int[] array;
		
		public IntArray(int[] array) {
			this.array = array;
		}

		@Override
		public int hashCode() {
			return Arrays.hashCode(array);
		}

		@Override
		public boolean equals(Object obj) {
			if (!IntArray.class.isInstance(obj)) {
				return false;
			}
			IntArray other = (IntArray) obj;
			return Arrays.equals(array, other.array);
		}
	}
	
	/**
	 * Empty constructor
	 *
	 */
	public KernelEstimatorMultiVariate() {
	}

	/**
	 * Initialise the estimator before passing any observations in.
	 *
	 * @param kernelWidth
	 * @param dimensions
	 */
	public void initialise(int dimensions, double kernelWidth) {
		this.dimensions = dimensions;
		this.suppliedKernelWidths = new double[dimensions];
		for (int d = 0; d < dimensions; d++) {
			this.suppliedKernelWidths[d] = kernelWidth;
		}
		finishInitialisation();
	}

	/**
	 * Initialise the estimator before passing any observations in.
	 *
	 * @param kernelWidths
	 */
	public void initialise(double[] kernelWidths) {
		dimensions = kernelWidths.length;
		this.suppliedKernelWidths = new double[dimensions];
		for (int d = 0; d < dimensions; d++) {
			this.suppliedKernelWidths[d] = kernelWidths[d];
		}
		finishInitialisation();
	}
	
	/**
	 * Private method to provide common variable initialisation
	 *
	 */
	private void finishInitialisation() {
		observations = null;
		observationsByHash = null;
		bins = null;
		useBins = true;
		usingIntegerIndexBins = true;
		mins = null;
		multipliers = null;
		rawData = null;
		totalObservations = 0;
		kernelWidthsInUse = new double[dimensions];
	}

	/**
	 * Each row of the data is an observation; each column of
	 *  the row is a new variable in the multivariate observation.
	 * 
	 * @param data
	 */
	@SuppressWarnings("unchecked")
	public void setObservations(double[][] data) {
		totalObservations = data.length;
		bins = new int[dimensions];
		mins = new double[dimensions];
		multipliers = new int[dimensions];
		int multiDimBins = 1;
		// If we are forcing n^2 comparisons, don't even consider using integer index bins.
		// Otherwise, we will assume we're doing this unless we find there would be too 
		//  many bins:
		usingIntegerIndexBins = !forceCompareToAll;
		useBins = true;
		// Work out the max and min for each bin
		for (int d = 0; d < dimensions; d++) {
			mins[d] = MatrixUtils.min(data, d);
			double max = MatrixUtils.max(data, d);
			double std = 0.0;
			if (normalise) {
				// Compute what the epsilonInUse should be here:
				//  it should expand with the standard deviation.
				// This saves us from normalising all of the incoming data points!
				std = MatrixUtils.stdDev(data, d);
				kernelWidthsInUse[d] = suppliedKernelWidths[d] * std;
			} else {
				kernelWidthsInUse[d] = suppliedKernelWidths[d];
			}
			bins[d] = (int) Math.ceil((max - mins[d]) / kernelWidthsInUse[d]);
			if (bins[d] == 0) {
				// This means the min and max are exactly the same:
				//  for our purposes this is akin to requiring one bin here.
				// If we leave it as zero bins, this will screw up
				//  our calculation of the number of multidimensional bins
				//  and cause errors later.
				bins[d] = 1;
			}
			multipliers[d] = multiDimBins;
			multiDimBins *= bins[d];
			if (usingIntegerIndexBins && (multiDimBins > MAX_NUMBER_INTEGER_BINS)) {
				// the number of bins has breached our maximum, so we'll use
				//  a hash table instead of an integer array
				usingIntegerIndexBins = false;
			}
			if (debug) {
				System.out.println("Dim: " + d + " => Max: " + max + ", min: " + mins[d] +
					", bins: " + bins[d] +
					(normalise ? ", std: " + std : "") +
					", eps: " + kernelWidthsInUse[d]);
			}
		}
		
		observations = null;
		observationsByHash = null;
		rawData = data; // Store the raw data regardless of our
						//  indexing approach - we may need to 
						//  access it for dynamic correlation exclusion
		if (usingIntegerIndexBins) {
			// multiDimBins is the number of multi-dimensional bins we have

			if (debug) {
				System.out.println("Multidimensional bins: " + multiDimBins);
			}
			
			// Create the bins
			observations = (Vector<TimeStampedObservation>[]) new Vector[multiDimBins];
			for (int v = 0; v < multiDimBins; v++) {
				observations[v] = new Vector<TimeStampedObservation>();
			}
			
			// Add each observation
			for (int t = 0; t < data.length; t++) {
				int bin = getMultiDimBinIndex(data[t]);
				// Create the time-stamped observation to store:
				TimeStampedObservation tso = new TimeStampedObservation(t, data[t]);
				observations[bin].add(tso);
			}

		} else {
			// Could be using either hash table style bins or no bins at all.
			// Work out how many bins we might be trawling through if we used bins:
			//  (need to take 3 ^ d, since we look at one bin either side for each dimension)
			int possibleBinsToTrawl = 1;
			for (int d = 0; d < dimensions; d++) {
				if (possibleBinsToTrawl > (Integer.MAX_VALUE / 3)) {
					// Next multiplication will overflow
					possibleBinsToTrawl = Integer.MAX_VALUE;
					break;
				}
				possibleBinsToTrawl *= 3;
			}
			
			// We might be forcing n^2 comparisons between all points, or 
			//  as a rule of thumb, check if there are less observations than possible bins to trawl
			if (forceCompareToAll || (totalObservations < possibleBinsToTrawl)) {
				useBins = false;
				
				// Simply store the raw data:
				// rawData = data; - this is already done above
				
			} else {
				// Using hash table
				useBins = true;
				
				if (debug) {
					System.out.println("Using array hash index bins");
				}
				
				// Create the hashtable
				observationsByHash = new Hashtable<IntArray, Vector<TimeStampedObservation>>();
				
				// Add each observation
				for (int t = 0; t < data.length; t++) {
					int[] multiDimBin = getMultiDimBinArray(data[t]);
					IntArray intArrayMultiDimBin = new IntArray(multiDimBin);
					Vector<TimeStampedObservation> hashedVector =
						(Vector<TimeStampedObservation>) observationsByHash.get(intArrayMultiDimBin);
					if (hashedVector == null) {
						hashedVector = new Vector<TimeStampedObservation>();
					}
					// Create the time-stamped observation to store:
					TimeStampedObservation tso = new TimeStampedObservation(t, data[t]);
					hashedVector.add(tso);
					// And put the observations for this multidimensional bin back in
					observationsByHash.put(intArrayMultiDimBin, hashedVector);
				}
			}
		}
		
		// Not sorting the observations, as this can only
		//  be practically done for a single variable.
	}
	
	/**
	 * Each row of the data is an observation; each column of
	 *  the row is a new variable in the multivariate observation.
	 * This method signature allows the user to call setObservations for
	 *  joint time series without combining them into a single joint time
	 *  series (we do the combining for them).
	 * 
	 * @param data1
	 * @param data2
	 * @throws Exception When the length of the two arrays of observations do not match.
	 */
	public void setObservations(double[][] data1, double[][] data2) throws Exception {
		int timeSteps = data1.length;
		if ((data1 == null) || (data2 == null)) {
			throw new Exception("Cannot have null data arguments");
		}
		if (data1.length != data2.length) {
			throw new Exception("Length of data1 (" + data1.length + ") is not equal to the length of data2 (" +
					data2.length + ")");
		}
		int data1Variables = data1[0].length;
		int data2Variables = data2[0].length;
		double[][] data = new double[timeSteps][data1Variables + data2Variables];
		for (int t = 0; t < timeSteps; t++) {
			System.arraycopy(data1[t], 0, data[t], 0, data1Variables);
			System.arraycopy(data2[t], 0, data[t], data1Variables, data2Variables);
		}
		// Now defer to the normal setObservations method
		setObservations(data);
	}
	
	/**
	 * This method signature allows the user to call setObservations for
	 *  joint time series without combining them into a single joint time
	 *  series (we do the combining for them).
	 * 
	 * @param data1
	 * @param data2
	 * @throws Exception When the length of the two arrays of observations do not match.
	 */
	public void setObservations(double[] data1, double[] data2) throws Exception {
		int timeSteps = data1.length;
		if ((data1 == null) || (data2 == null)) {
			throw new Exception("Cannot have null data arguments");
		}
		if (data1.length != data2.length) {
			throw new Exception("Length of data1 (" + data1.length + ") is not equal to the length of data2 (" +
					data2.length + ")");
		}
		double[][] data = new double[timeSteps][2];
		MatrixUtils.copyIntoColumn(data, 0, data1);
		MatrixUtils.copyIntoColumn(data, 1, data2);
		// Now defer to the normal setObservations method
		setObservations(data);
	}
	
	/**
	 * Return the kernel estimate for the probability without any dynamic
	 *  correlation exlcusion
	 * 
	 * @param observation
	 * @return
	 */
	public double getProbability(double[] observation) {
		if (useBins) {
			return getProbabilityFromBins(observation, 0, false);
		} else {
			return getProbabilityComparingToAll(observation, 0, false);
		}
	}

	/**
	 * Return the kernel estimate for the probability with dynamic
	 *  correlation exlcusion if it has been swithced on
	 * 
	 * @param observation
	 * @param timeStep
	 * @return
	 */
	public double getProbability(double[] observation, int timeStep) {
		if (useBins) {
			return getProbabilityFromBins(observation, timeStep, excludeDynamicCorrelations);
		} else {
			return getProbabilityComparingToAll(observation, timeStep, excludeDynamicCorrelations);
		}
	}

	/**
	 * Return the kernel estimate for the joint probability without any dynamic
	 *  correlation exlcusion
	 *  
	 * @param observation1
	 * @param observation2
	 * @return
	 */
	public double getProbability(double[] observation1, double[] observation2) {
		// Make a joint array
		double[] jointObservation = new double[observation1.length + observation2.length];
		System.arraycopy(observation1, 0, jointObservation, 0, observation1.length);
		System.arraycopy(observation2, 0, jointObservation, observation1.length, observation2.length);
		return getProbability(jointObservation);
	}

	/**
	 * Return the kernel estimate for the joint probability with dynamic
	 *  correlation exlcusion if it has been swithced on
	 *  
	 * @param observation1
	 * @param observation2
	 * @param timeStep
	 * @return
	 */
	public double getProbability(double[] observation1, double[] observation2, int timeStep) {
		// Make a joint array
		double[] jointObservation = new double[observation1.length + observation2.length];
		System.arraycopy(observation1, 0, jointObservation, 0, observation1.length);
		System.arraycopy(observation2, 0, jointObservation, observation1.length, observation2.length);
		return getProbability(jointObservation, timeStep);
	}

	/**
	 * Return the kernel count (i think this is the correlation integral effectively?)
	 *  for this vector without any dynamic correlation exlcusion
	 * 
	 * @param observation
	 * @return
	 */
	public int getCount(double[] observation) {
		KernelCount kernelCount;
		if (useBins) {
			kernelCount = getCountFromBins(observation, 0, false);
		} else {
			kernelCount = getCountComparingToAll(observation, 0, false);
		}
		return kernelCount.count;
	}

	/**
	 * <p>Return the kernel count for the probability with dynamic
	 *  correlation exlcusion if it has been switched on.</p>
	 *  
	 * <p>The user should be aware that a call to getProbability(double[], int)
	 * will not return the same value as getCount(double[], int) / totalObservations,
	 * since dynamic correlation exclusion could be swithced on. To achieve this, the
	 * user could call getCompleteKernelCount(double[], int[], boolean) and use the
	 * returned count and totalObservationsForCount in the returned KernelCount object.
	 * </p>
	 * 
	 * @param observation
	 * @param timeStep
	 * @return
	 */
	public int getCount(double[] observation, int timeStep) {
		KernelCount kernelCount;
		if (useBins) {
			kernelCount = getCountFromBins(observation, timeStep, excludeDynamicCorrelations);
		} else {
			kernelCount = getCountComparingToAll(observation, timeStep, excludeDynamicCorrelations);
		}
		return kernelCount.count;
	}

	/**
	 * Return the kernel count for the probability with dynamic
	 *  correlation exlcusion if it has been switched on
	 * 
	 * @param observation
	 * @param timeStep
	 * @param giveListOfCorrelatedPoints
	 * @return
	 */
	public KernelCount getCompleteKernelCount(double[] observation, int timeStep,
			boolean giveListOfCorrelatedPoints) {
		KernelCount kernelCount;
		if (useBins) {
			kernelCount = getCountFromBins(observation, timeStep,
					excludeDynamicCorrelations, giveListOfCorrelatedPoints);
		} else {
			kernelCount = getCountComparingToAll(observation, timeStep,
					excludeDynamicCorrelations, giveListOfCorrelatedPoints);
		}
		return kernelCount;
	}

	/**
	 * Return the kernel count for the joint vector without any dynamic
	 *  correlation exlcusion
	 *  
	 * @param observation1
	 * @param observation2
	 * @return
	 */
	public int getCount(double[] observation1, double[] observation2) {
		// Make a joint array
		double[] jointObservation = new double[observation1.length + observation2.length];
		System.arraycopy(observation1, 0, jointObservation, 0, observation1.length);
		System.arraycopy(observation2, 0, jointObservation, observation1.length, observation2.length);
		return getCount(jointObservation);
	}

	/**
	 * Return the kernel estimate for the joint probability with dynamic
	 *  correlation exlcusion if it has been swithced on
	 *  
	 * @param observation1
	 * @param observation2
	 * @param timeStep
	 * @return
	 */
	public int getCount(double[] observation1, double[] observation2, int timeStep) {
		// Make a joint array
		double[] jointObservation = new double[observation1.length + observation2.length];
		System.arraycopy(observation1, 0, jointObservation, 0, observation1.length);
		System.arraycopy(observation2, 0, jointObservation, observation1.length, observation2.length);
		return getCount(jointObservation, timeStep);
	}

	/**
	 * Return the kernel count object for the joint probability with dynamic
	 *  correlation exlcusion if it has been swithced on
	 *  
	 * @param observation1
	 * @param observation2
	 * @param timeStep
	 * @return
	 */
	public KernelCount getCompleteKernelCount(double[] observation1, double[] observation2,
			int timeStep, boolean giveListOfCorrelatedPoints) {
		// Make a joint array
		double[] jointObservation = new double[observation1.length + observation2.length];
		System.arraycopy(observation1, 0, jointObservation, 0, observation1.length);
		System.arraycopy(observation2, 0, jointObservation, observation1.length, observation2.length);
		return getCompleteKernelCount(jointObservation, timeStep, giveListOfCorrelatedPoints);
	}

	/**
	 * Give kernel estimated count for this observation.
	 *  Do this using bins, so we only need to compare to neighbouring bins 
	 *   and of course add in the observations in this bin.
	 *  Achieves this by trawling through neighbouring bins in a recursive fashion.
	 * 
	 * @param observation
	 * @param timeStep of the given observation (required for dynamic correlation exclusion)
	 * @param dynCorrExclusion whether to use dynamic correlation exclusion
	 * @return KernelCount object with the count and total observations the count was taken from
	 */
	private KernelCount getCountFromBins(double[] observation, int timeStep, 
			boolean dynCorrExclusion) {
		return getCountFromBins(observation, timeStep, dynCorrExclusion, false);
	}

	/**
	 * Give kernel estimated count for this observation.
	 *  Do this using bins, so we only need to compare to neighbouring bins 
	 *   and of course add in the observations in this bin.
	 *  Achieves this by trawling through neighbouring bins in a recursive fashion.
	 * 
	 * @param observation
	 * @param timeStep of the given observation (required for dynamic correlation exclusion)
	 * @param dynCorrExclusion whether to use dynamic correlation exclusion
	 * @param giveListOfCorrelatedPoints include the list of time indices whose vectors
	 *  were correlated with the given vector
	 * @return KernelCount object with the count and total observations the count was taken from
	 */
	protected KernelCount getCountFromBins(double[] observation, int timeStep, 
			boolean dynCorrExclusion, boolean giveListOfCorrelatedPoints) {
		// First count the number of observations in the same bin
		int count;
		boolean[] isCorrelatedWithArgument = null;
		if (giveListOfCorrelatedPoints) {
			isCorrelatedWithArgument = new boolean[totalObservations];
		}
		int multiDimBin = 0;
		int[] multiDimBinArray = null;
		if (usingIntegerIndexBins) {
			multiDimBin = getMultiDimBinIndex(observation);
			count = observations[multiDimBin].size();
			if (giveListOfCorrelatedPoints || makeCorrelatedPointAddedCallback) {
				for (int i = 0; i < observations[multiDimBin].size(); i++) {
					TimeStampedObservation tso = (TimeStampedObservation) observations[multiDimBin].elementAt(i);
					if (giveListOfCorrelatedPoints) {
						isCorrelatedWithArgument[tso.timeStep] = true;
					}
					if (makeCorrelatedPointAddedCallback) {
						correlatedPointAddedCallback(tso.timeStep);
					}
				}
			}
		} else {
			multiDimBinArray = getMultiDimBinArray(observation);
			IntArray intArrayMultiDimBin = new IntArray(multiDimBinArray);
			Vector<TimeStampedObservation> observationsInThisBin =
				(Vector<TimeStampedObservation>) observationsByHash.get(intArrayMultiDimBin);
			// Shouldn't need to check (observationsInThisBin != null) here since 
			// at least this observation itself should have been added in here at some point.
			// Will need to change this though if we ever start using the kernel estimator to check
			// probabilities for new vectors that were not previously added to the observations.
			count = observationsInThisBin.size();
			if (giveListOfCorrelatedPoints || makeCorrelatedPointAddedCallback) {
				for (int i = 0; i < observationsInThisBin.size(); i++) {
					TimeStampedObservation tso = (TimeStampedObservation) observationsInThisBin.elementAt(i);
					if (giveListOfCorrelatedPoints) {
						isCorrelatedWithArgument[tso.timeStep] = true;
					}
					if (makeCorrelatedPointAddedCallback) {
						correlatedPointAddedCallback(tso.timeStep);
					}
				}
			}
		}
		
		int totalTimePointsCompared = totalObservations;
		
		if (dynCorrExclusion) {
			// Need to remove any observations that were *closer* than timeProximityForDynamicCorrelationExclusion
			int closeTimePointsToCompare = (timeStep >= timeProximityForDynamicCorrelationExclusion) ?
					timeProximityForDynamicCorrelationExclusion - 1: timeStep;
			closeTimePointsToCompare += (totalObservations - timeStep >= timeProximityForDynamicCorrelationExclusion) ?
					timeProximityForDynamicCorrelationExclusion - 1: totalObservations - timeStep - 1;
			closeTimePointsToCompare++; // Add one for comparison to self
			totalTimePointsCompared -= closeTimePointsToCompare;
			// Use a heuristic to select best way to eliminate dynamic correlations here:
			int countToRemove = 0;
			if (closeTimePointsToCompare * dimensions < count) {
				// Check all the neighbouring points in time to see if
				//  they were in the same bin and have been counted
				for (int t = Math.max(0,timeStep-timeProximityForDynamicCorrelationExclusion+1);
						t < Math.min(totalObservations,timeStep+timeProximityForDynamicCorrelationExclusion); t++) {
					int otherMultiDimBin = getMultiDimBinIndex(rawData[t]);
					if (otherMultiDimBin == multiDimBin) {
						countToRemove++;
						if (giveListOfCorrelatedPoints) {
							isCorrelatedWithArgument[t] = false;
						}
						if (makeCorrelatedPointAddedCallback) {
							correlatedPointRemovedCallback(t);
						}
					}
				}
			} else {
				// check all the points in this bin 
				if (usingIntegerIndexBins) {
					for (int i = 0; i < observations[multiDimBin].size(); i++) {
						TimeStampedObservation tso = (TimeStampedObservation) observations[multiDimBin].elementAt(i);
						if (Math.abs(tso.timeStep - timeStep) < timeProximityForDynamicCorrelationExclusion) {
							countToRemove++;
							if (giveListOfCorrelatedPoints) {
								isCorrelatedWithArgument[tso.timeStep] = false;
							}
							if (makeCorrelatedPointAddedCallback) {
								correlatedPointRemovedCallback(tso.timeStep);
							}
						}
					}
				} else {
					Vector<TimeStampedObservation> observationsInThisBin =
						(Vector<TimeStampedObservation>) observationsByHash.get(multiDimBinArray);
					if (observationsInThisBin != null) {
						for (int i = 0; i < observationsInThisBin.size(); i++) {
							TimeStampedObservation tso = (TimeStampedObservation) observationsInThisBin.elementAt(i);
							if (Math.abs(tso.timeStep - timeStep) < timeProximityForDynamicCorrelationExclusion) {
								countToRemove++;
								if (giveListOfCorrelatedPoints) {
									isCorrelatedWithArgument[tso.timeStep] = false;
								}
								if (makeCorrelatedPointAddedCallback) {
									correlatedPointRemovedCallback(tso.timeStep);
								}
							}
						}
					}
				}
			}
			count -= countToRemove;
		}

		// Now look for correlated points in neighbouring bins
		int[] currentNeighbourBinIndices = new int[dimensions];
		// Run for one bin lower (assuming it exists)
		count += addCount(observation, currentNeighbourBinIndices, 0,
				getBinIndex(observation[0], 0) - 1, false,
				dynCorrExclusion, timeStep, isCorrelatedWithArgument);
		// Run along this bin for the moment, changing one of the other indices
		count += addCount(observation,currentNeighbourBinIndices, 0,
				getBinIndex(observation[0], 0), true,
				dynCorrExclusion, timeStep, isCorrelatedWithArgument);
		// Run for the next bin up (assuming it exists)
		count += addCount(observation, currentNeighbourBinIndices, 0,
				getBinIndex(observation[0], 0) + 1, false,
				dynCorrExclusion, timeStep, isCorrelatedWithArgument);
		
		KernelCount kernelCount = new KernelCount(count, totalTimePointsCompared, isCorrelatedWithArgument);
		return kernelCount;
	}

	/**
	 * Give kernel estimated probability for this observation.
	 *  Do this using bins, so we only need to compare to neighbouring bins 
	 *   and of course add in the observations in this bin.
	 *  Achieves this by trawling through neighbouring bins in a recursive fashion.
	 * 
	 * @param observation
	 * @param timeStep of the given observation (required for dynamic correlation exclusion)
	 * @param dynCorrExclusion whether to use dynamic correlation exclusion
	 * @return
	 */
	private double getProbabilityFromBins(double[] observation, int timeStep, 
			boolean dynCorrExclusion) {
		
		KernelCount kernelCount = getCountFromBins(observation, timeStep, dynCorrExclusion);
		
		return (double) kernelCount.count / (double) kernelCount.totalObservationsForCount;
	}
	
	/**
	 * Give kernel estimated count for this observation, by comparing to all
	 *  given observations.
	 * 
	 * @param observation
	 * @param timeStep of the given observation (required for dynamic correlation exclusion)
	 * @param dynCorrExclusion whether to use dynamic correlation exclusion
	 * @return KernelCount object with the count and total observations the count was taken from
	 */
	private KernelCount getCountComparingToAll(double[] observation, int timeStep, 
			boolean dynCorrExclusion) {
		return getCountComparingToAll(observation, timeStep, dynCorrExclusion, false);
	}


	/**
	 * Give kernel estimated count for this observation, by comparing to all
	 *  given observations.
	 * 
	 * @param observation
	 * @param timeStep of the given observation (required for dynamic correlation exclusion)
	 * @param dynCorrExclusion whether to use dynamic correlation exclusion
	 * @return KernelCount object with the count and total observations the count was taken from
	 */
	protected KernelCount getCountComparingToAll(double[] observation, int timeStep, 
			boolean dynCorrExclusion, boolean giveListOfCorrelatedPoints) {
		int count = 0;
		boolean[] isCorrelatedWithArgument = null;
		if (giveListOfCorrelatedPoints) {
			isCorrelatedWithArgument = new boolean[totalObservations];
		}
		for (int t = 0; t < totalObservations; t++) {
			if (!dynCorrExclusion || 
					(Math.abs(t - timeStep) >= timeProximityForDynamicCorrelationExclusion)) {
				// Only add in the contribution from this point
				//  when we're not doing dynamic correlation exclusion
				//  or if it's outside the exclusion zone
				int thisStepKernelValue = stepKernel(observation, rawData[t]);
				count += thisStepKernelValue;
				if (giveListOfCorrelatedPoints) {
					isCorrelatedWithArgument[t] = (thisStepKernelValue > 0);
				}
				if (makeCorrelatedPointAddedCallback && (thisStepKernelValue > 0)) {
					correlatedPointAddedCallback(t);
				}
			}
		}
		int totalTimePointsCompared = totalObservations;
		if (dynCorrExclusion) {
			// Needed to remove any observations that were *closer* than timeProximityForDynamicCorrelationExclusion
			//  from the total of time points that were compared.
			int closeTimePoints = (timeStep >= timeProximityForDynamicCorrelationExclusion) ?
					timeProximityForDynamicCorrelationExclusion - 1: timeStep;
			closeTimePoints += (totalObservations - timeStep >= timeProximityForDynamicCorrelationExclusion) ?
					timeProximityForDynamicCorrelationExclusion - 1: totalObservations - timeStep - 1;
			closeTimePoints++; // Add one for comparison to self
			totalTimePointsCompared -= closeTimePoints;
		}
		KernelCount kernelCount = new KernelCount(count, totalTimePointsCompared, isCorrelatedWithArgument);
		return kernelCount;
	}

	/**
	 * Give kernel estimated probability for this observation, by comparing to all
	 *  given observations.
	 * 
	 * @param observation
	 * @param timeStep of the given observation (required for dynamic correlation exclusion)
	 * @param dynCorrExclusion whether to use dynamic correlation exclusion
	 * @return
	 */
	private double getProbabilityComparingToAll(double[] observation, int timeStep, 
			boolean dynCorrExclusion) {
		
		KernelCount kernelCount = getCountComparingToAll(observation, timeStep, dynCorrExclusion);
		return (double) kernelCount.count / (double) kernelCount.totalObservationsForCount;
	}
	
	private int getBinIndex(double value, int dimension) {
		int bin = (int) Math.floor((value - mins[dimension]) / kernelWidthsInUse[dimension]);
		// Check for any rounding errors on the bin assignment:
		if (bin >= bins[dimension]) {
			bin = bins[dimension] - 1;
		}
		if (bin < 0) {
			bin = 0;
		}
		return bin;
	}
	
	public int getMultiDimBinIndex(double[] value) {
		int multiDimBin = 0;
		for (int d = 0; d < dimensions; d++) {
			multiDimBin += getBinIndex(value[d], d) * multipliers[d];
		}
		return multiDimBin;
	}

	public int[] getMultiDimBinArray(double[] value) {
		int[] multiDimBinArray = new int[dimensions];
		for (int d = 0; d < dimensions; d++) {
			multiDimBinArray[d] = getBinIndex(value[d], d);
		}
		return multiDimBinArray;
	}

	public int getMultiDimBinIndexFromSingles(int[] singleBinIndices) {
		int multiDimBin = 0;
		for (int d = 0; d < dimensions; d++) {
			multiDimBin += singleBinIndices[d] * multipliers[d];
		}
		return multiDimBin;
	}

	/**
	 * Count the number of data points within epsilon, 
	 * for multi dimensional bin indices as specified
	 * 
	 * @param neighbourBinIndices
	 * @param currentIndex
	 * @param currentIndexValue
	 * @param possibleCentralBin
	 * @param dynCorrExclusion
	 * @param timeStep
	 * @param isCorrelatedWithArgument set if the vector at each given time point is included 
	 * @return
	 */
	private int addCount(double[] observation, int[] neighbourBinIndices,
			int currentIndex, int currentIndexValue, boolean possibleCentralBin,
			boolean dynCorrExclusion, int timeStep, boolean[] isCorrelatedWithArgument) {
		
		int count = 0;
		
		// Check that the current index is in range
		if (currentIndexValue < 0) {
			return 0;
		}
		if (currentIndexValue >= bins[currentIndex]) {
			return 0;
		}
		neighbourBinIndices[currentIndex] = currentIndexValue;
		
		if (currentIndex == dimensions - 1) {
			// this is the last dimension, need to check the neighbouring bin here
			if (possibleCentralBin) {
				// This is definitely the central bin, don't
				//  count it, we've already done that outside
				return 0;
			}
			if (usingIntegerIndexBins) {
				int multiDimIndex = getMultiDimBinIndexFromSingles(neighbourBinIndices);
				for (int i = 0; i < observations[multiDimIndex].size(); i++) {
					TimeStampedObservation tso = (TimeStampedObservation) observations[multiDimIndex].elementAt(i);
					if (!dynCorrExclusion || 
							(Math.abs(tso.timeStep - timeStep) >= timeProximityForDynamicCorrelationExclusion)) {
						// Only add in the contribution from this point
						//  when we're not doing dynamic correlation exclusion
						//  or if it's outside the exclusion zone
						int thisStepKernelValue = stepKernel(observation, tso.observation);
						count += thisStepKernelValue;
						if (isCorrelatedWithArgument != null) {
							isCorrelatedWithArgument[tso.timeStep] = (thisStepKernelValue > 0);
						}
						if (makeCorrelatedPointAddedCallback && (thisStepKernelValue > 0)) {
							correlatedPointAddedCallback(tso.timeStep);
						}
					}
				}
			} else {
				IntArray intArrayBinNeighbourIndices = new IntArray(neighbourBinIndices);
				Vector<TimeStampedObservation> observationsInThisBin =
					(Vector<TimeStampedObservation>) observationsByHash.get(intArrayBinNeighbourIndices);
				if (observationsInThisBin != null) {
					for (int i = 0; i < observationsInThisBin.size(); i++) {
						TimeStampedObservation tso = (TimeStampedObservation) observationsInThisBin.elementAt(i);
						if (!dynCorrExclusion || 
								(Math.abs(tso.timeStep - timeStep) >= timeProximityForDynamicCorrelationExclusion)) {
							// Only add in the contribution from this point
							//  when we're not doing dynamic correlation exclusion
							//  or if it's outside the exclusion zone
							int thisStepKernelValue = stepKernel(observation, tso.observation);
							count += thisStepKernelValue;
							if (isCorrelatedWithArgument != null) {
								isCorrelatedWithArgument[tso.timeStep] = (thisStepKernelValue > 0);
							}
							if (makeCorrelatedPointAddedCallback && (thisStepKernelValue > 0)) {
								correlatedPointAddedCallback(tso.timeStep);
							}
						}
					}
				}
			}
		} else {
			// pass on the responsibility to do the checks:
			// Run for one index lower (assuming it exists)
			count += addCount(observation, neighbourBinIndices, currentIndex+1,
					getBinIndex(observation[currentIndex+1], currentIndex+1) - 1, false,
					dynCorrExclusion, timeStep, isCorrelatedWithArgument);
			// Run for this index
			count += addCount(observation,neighbourBinIndices, currentIndex+1,
					getBinIndex(observation[currentIndex+1], currentIndex+1), possibleCentralBin,
					dynCorrExclusion, timeStep, isCorrelatedWithArgument);
			// Run for the next index up (assuming it exists)
			count += addCount(observation, neighbourBinIndices, currentIndex+1,
					getBinIndex(observation[currentIndex+1], currentIndex+1) + 1, false,
					dynCorrExclusion, timeStep, isCorrelatedWithArgument);
		}
		
		return count;
	}
	
	/**
	 * Return the step kernel for the two data points
	 * 
	 * @param observation
	 * @param candidate
	 * @return
	 */
	public int stepKernel(double[] observation, double[] candidate) {
		for (int d = 0; d < dimensions; d++) {
			if (Math.abs(observation[d] - candidate[d]) > kernelWidthsInUse[d]) {
				return 0;
			}
		}
		// The candidate is within epsilon of the observation
		return 1;
	}

	/**
	 * Return the kernel widths that are actually in use in the
	 * calculator (i.e. those after any normalisation is applied)
	 * 
	 * @return the kernelWidthsInUse
	 */
	public double[] getKernelWidthsInUse() {
		return kernelWidthsInUse;
	}

	public void setNormalise(boolean normalise) {
		this.normalise = normalise;
	}
	
	public boolean isNormalise() {
		return normalise;
	}

	public void setDynamicCorrelationExclusion(int timeWindow) {
		excludeDynamicCorrelations = true;
		timeProximityForDynamicCorrelationExclusion = timeWindow;
	}
	
	public void clearDynamicCorrelationExclusion() {
		excludeDynamicCorrelations = false;
	}
	
	public boolean isExcludeDynamicCorrelations() {
		return excludeDynamicCorrelations;
	}

	public boolean isForceCompareToAll() {
		return forceCompareToAll;
	}

	/**
	 * Only recommended to set this for debugging purposes
	 * 
	 * @param forceCompareToAll
	 */
	public void setForceCompareToAll(boolean forceCompareToAll) {
		this.forceCompareToAll = forceCompareToAll;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}
	
	/**
	 * A callback for where a correlated point is found.
	 * Used in child classes, but not here since we have set makeCorrelatedPointAddedCallback to false.
	 *
	 */
	protected void correlatedPointAddedCallback(int correlatedTimeStep) {
		
	}

	/**
	 * A callback for where a correlated point is removed due to dynamic correlated exclusion.
	 * Used in child classes, but not here since we have set makeCorrelatedPointAddedCallback to false.
	 *
	 */
	protected void correlatedPointRemovedCallback(int removedCorrelatedTimeStep) {
		
	}

	/**
	 * Provide an implementation of the clone() method.
	 * This does not deeply copy all of the underlying data, just providing
	 *  a copy of the references to it all.
	 * This is enough to protect the integrity of the kernel estimator
	 *  however if the clone is supplied different data (though the 
	 *  clone should not alter the data).
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	protected Object clone() throws CloneNotSupportedException {
		return super.clone();
	}
}
