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

package infodynamics.utils;

import java.util.Collection;
import java.util.PriorityQueue;
import java.util.Vector;


/**
 * Class for fast neighbour searching
 *  in a <b>single dimensional</b> variable.
 * Instantiates a sorted array for this purpose.
 * Norms for the nearest neighbour searches are the max norm or Euclidean
 *  norm (squared) within each variable.
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class UnivariateNearestNeighbourSearcher extends NearestNeighbourSearcher {

	/**
	 * Cached reference to the original data set
	 */
	protected double[] originalDataSet;
	/**
	 * Number of samples in the data
	 */
	protected int numObservations = 0;
	/**
	 * An array of indices to the data in originalDataSet,
	 *  sorted in order (min to max).
	 */
	protected int[] sortedArrayIndices = null;
	
	/**
	 * An array of indices of where each data point
	 *  in originalDataSet lies in the sorted array
	 */
	protected int[] indicesInSortedArray = null;
	
	public UnivariateNearestNeighbourSearcher(double[][] data) throws Exception {
		// Ideally we would not call the constructor until after the following check,
		//  but the constructor must come first in Java.
		this(MatrixUtils.selectColumn(data, 0));
		if (data[0].length != 1) {
			throw new Exception("Cannot define UnivariateNearestNeighbourSearcher for multivariate data");
		}
	}
	
	public UnivariateNearestNeighbourSearcher(double[] data) throws Exception {
		this.originalDataSet = data;
		numObservations = data.length;
		if (numObservations <= 1) {
			throw new Exception("Nearest neighbour search is poorly defined for <=1 data point");
		}
		
		// Sort the original data sets in each dimension:
		double[][] dataWithIndices = new double[numObservations][2];
		// Record original time indices:
		for (int t = 0; t < numObservations; t++) {
			dataWithIndices[t][0] = data[t];
			dataWithIndices[t][1] = t;
		}
		// Sort the data:
		java.util.Arrays.sort(dataWithIndices, FirstIndexComparatorDouble.getInstance());
		// And extract the sorted indices, and references
		//  to where each original sample lies in the sorted array
		sortedArrayIndices = new int[numObservations];
		indicesInSortedArray = new int[numObservations];
		for (int t = 0; t < numObservations; t++) {
			sortedArrayIndices[t] = (int) dataWithIndices[t][1];
			indicesInSortedArray[(int) dataWithIndices[t][1]] = t;
		}
	}

	/**
	 * Computed the configured norm between the unidimensional variables.
	 * Hoping this method is inlined by the JVM, but haven't checked this.
	 * 
	 * @param x1 data point 1
	 * @param x2 data point 2
	 * @return the norm
	 */
	protected static double norm(double x1, double x2, int normTypeToUse) {
		switch (normTypeToUse) {
		case EuclideanUtils.NORM_MAX_NORM:
			return Math.abs(x1-x2);
		// case EuclideanUtils.NORM_EUCLIDEAN_SQUARED:
		default:
			double difference = x1 - x2;
			return difference * difference;
		}
	}
	
	/**
	 * Return the node which is the nearest neighbour for a given
	 *  sample index in the data set. The node itself is 
	 *  excluded from the search.
	 * 
	 */
	public NeighbourNodeData findNearestNeighbour(int sampleIndex) {
		// Find where this node sits in the sorted array:
		int indexInSortedArray = indicesInSortedArray[sampleIndex];
		if (indexInSortedArray == 0) {
			// There is only one candidate for nearest neighbour
			// Assumes we have more than 1 data point -- this is
			//  checked in the constructor for us.
			double theNorm = norm(originalDataSet[sampleIndex],
						originalDataSet[sortedArrayIndices[1]], normTypeToUse);
			return new NeighbourNodeData(sortedArrayIndices[1],
					new double[] {theNorm}, theNorm);
		} else if (indexInSortedArray == numObservations - 1) {
			// There is only one candidate for nearest neighbour
			// Assumes we have more than 1 data point -- this is
			//  checked in the constructor for us.
			double theNorm = norm(originalDataSet[sampleIndex],
						originalDataSet[sortedArrayIndices[numObservations - 2]],
						normTypeToUse);
			return new NeighbourNodeData(sortedArrayIndices[numObservations - 2],
					new double[] {theNorm}, theNorm);
		} else {
			// We need to check candidates on both sides of the data point:
			double normAbove = norm(originalDataSet[sampleIndex],
					originalDataSet[sortedArrayIndices[indexInSortedArray+1]],
					normTypeToUse);
			double normBelow = norm(originalDataSet[sampleIndex],
					originalDataSet[sortedArrayIndices[indexInSortedArray-1]],
					normTypeToUse);
			if (normAbove < normBelow) {
				return new NeighbourNodeData(sortedArrayIndices[indexInSortedArray+1],
						new double[] {normAbove}, normAbove);
			} else {
				return new NeighbourNodeData(sortedArrayIndices[indexInSortedArray-1],
						new double[] {normAbove}, normAbove);
			}
		}
	}
	
	/**
	 * Return the K nodes which are the K nearest neighbours for a given
	 *  sample index in the data set. The node itself is 
	 *  excluded from the search.
	 * Nearest neighbour function to compare to r is the specified norm.
	 */
	public PriorityQueue<NeighbourNodeData>
			findKNearestNeighbours(int K, int sampleIndex) throws Exception {
		
		if (numObservations <= K) {
			throw new Exception("Not enough data points for a K nearest neighbours search");
		}
		
		// Find where this node sits in the sorted array:
		int indexInSortedArray = indicesInSortedArray[sampleIndex];
		// Initialise the nearest neighbours above and below this data point,
		//  storing a -1 for the indices where there are none left on this side:
		int lowerCandidate = (indexInSortedArray == 0) ? -1 : indexInSortedArray - 1;
		int upperCandidate = (indexInSortedArray == numObservations - 1) ?
								-1 : indexInSortedArray + 1;
		
		PriorityQueue<NeighbourNodeData> pq = new PriorityQueue<NeighbourNodeData>(K);
		for (int k = 0; k < K; k++) {
			// Select the (k+1)th nearest neighbour
			// Check norms for candidates on both sides of the data point.
			//  (Their must be at least one valid candidate (i.e. not index -1)
			//   since we have previously checked there were at least K+1 data points)
			double normAbove =  (upperCandidate == -1) ?
					Double.POSITIVE_INFINITY :
						norm(originalDataSet[sampleIndex],
							originalDataSet[sortedArrayIndices[upperCandidate]],
							normTypeToUse);
			double normBelow = (lowerCandidate == -1) ?
					Double.POSITIVE_INFINITY :
						norm(originalDataSet[sampleIndex],
							originalDataSet[sortedArrayIndices[lowerCandidate]],
							normTypeToUse);
			NeighbourNodeData nextNearest;
			if (normAbove < normBelow) {
				nextNearest = new NeighbourNodeData(sortedArrayIndices[upperCandidate],
								new double[] {normAbove}, normAbove);
				// Advance the upper candidate
				upperCandidate = (upperCandidate == numObservations - 1) ?
									-1 : upperCandidate + 1;
			} else {
				nextNearest = new NeighbourNodeData(sortedArrayIndices[lowerCandidate],
								new double[] {normBelow}, normBelow);
				// Advance the lower candidate
				//lowerCandidate = (lowerCandidate == 0) ?
				//					-1 : lowerCandidate - 1;
				// Faster:
				lowerCandidate--;
			}
			// And add this next nearest neighbour to the PQ:
			pq.add(nextNearest);
		}
		
		return pq;
	}

	/**
	 * Return the K nodes which are the K nearest neighbours for a given
	 *  sample index in the data set. Nodes within dynCorrExclTime time
	 *  points are excluded from the search.
	 * Nearest neighbour function to compare to r is the specified norm.
	 */
	public PriorityQueue<NeighbourNodeData>
			findKNearestNeighbours(int K, int sampleIndex,
					int dynCorrExclTime) throws Exception {
		
		if (numObservations <= K + 2*dynCorrExclTime) {
			throw new Exception("Not enough data points for a K nearest neighbours search" +
					" with dynamic exclusion window of " + dynCorrExclTime + " points either side");
		}
		
		// Find where this node sits in the sorted array:
		int indexInSortedArray = indicesInSortedArray[sampleIndex];
		// Initialise the nearest neighbours above and below this data point,
		//  storing an extreme value for the indices where there are none left on this side:
		int lowerCandidate;
		int indexOfLowerCandidate = -1;
		for (lowerCandidate = indexInSortedArray - 1;
				lowerCandidate >= 0;
				lowerCandidate--) {
			indexOfLowerCandidate = sortedArrayIndices[lowerCandidate];
			if (Math.abs(sampleIndex - indexOfLowerCandidate) > dynCorrExclTime) {
				// This sample is outside the dynamic correlation exclusion window
				break;
			}
		}
		// Postcondition: lowerCandidate points to the nearest neighbour
		//  below this point, which is outside of the dynamic correlation exclusion
		//  window; otherwise lowerCandidate == -1.
		int upperCandidate;
		int indexOfUpperCandidate = -1;
		for (upperCandidate = indexInSortedArray + 1;
				upperCandidate <= numObservations - 1;
				upperCandidate++) {
			indexOfUpperCandidate = sortedArrayIndices[upperCandidate];
			if (Math.abs(sampleIndex - indexOfUpperCandidate) > dynCorrExclTime) {
				// This sample is outside the dynamic correlation exclusion window
				break;
			}
		}
		// Postcondition: upperCandidate points to the nearest neighbour
		//  above this point, which is outside of the dynamic correlation exclusion
		//  window; otherwise upperCandidate == numObservations.
		
		PriorityQueue<NeighbourNodeData> pq = new PriorityQueue<NeighbourNodeData>(K);
		for (int k = 0; k < K; k++) {
			// Select the (k+1)th nearest neighbour
			// Check norms for candidates on both sides of the data point.
			//  (Their must be at least one valid candidate
			//   since we have previously checked there were at least K+1 data points
			//   plus the dynamic exclusion window size)
			double normAbove =  (upperCandidate == numObservations) ?
					Double.POSITIVE_INFINITY :
						norm(originalDataSet[sampleIndex],
							originalDataSet[indexOfUpperCandidate],
							normTypeToUse);
			double normBelow = (lowerCandidate == -1) ?
					Double.POSITIVE_INFINITY :
						norm(originalDataSet[sampleIndex],
							originalDataSet[indexOfLowerCandidate],
							normTypeToUse);
			NeighbourNodeData nextNearest;
			if (normAbove < normBelow) {
				nextNearest = new NeighbourNodeData(indexOfUpperCandidate,
								new double[] {normAbove}, normAbove);
				// Advance the upper candidate
				for (upperCandidate++;
						upperCandidate <= numObservations - 1;
						upperCandidate++) {
					indexOfUpperCandidate = sortedArrayIndices[upperCandidate];
					if (Math.abs(sampleIndex - indexOfUpperCandidate) > dynCorrExclTime) {
						// This sample is outside the dynamic correlation exclusion window
						break;
					}
				}
				// Postcondition: upperCandidate points to the nearest neighbour
				//  above this point, which is outside of the dynamic correlation exclusion
				//  window; otherwise upperCandidate == numObservations.
			} else {
				nextNearest = new NeighbourNodeData(indexOfLowerCandidate,
								new double[] {normBelow}, normBelow);
				// Advance the lower candidate
				for (lowerCandidate--;
						lowerCandidate >= 0;
						lowerCandidate--) {
					indexOfLowerCandidate = sortedArrayIndices[lowerCandidate];
					if (Math.abs(sampleIndex - indexOfLowerCandidate) > dynCorrExclTime) {
						// This sample is outside the dynamic correlation exclusion window
						break;
					}
				}
				// Postcondition: lowerCandidate points to the nearest neighbour
				//  below this point, which is outside of the dynamic correlation exclusion
				//  window; otherwise lowerCandidate == -1.
			}
			// And add this next nearest neighbour to the PQ:
			pq.add(nextNearest);
		}
		
		return pq;
	}

	/**
	 * Count the number of points within norm r for a given
	 *  sample index in the data set. The node itself is 
	 *  excluded from the search.
	 * Nearest neighbour function to compare to r is the specified norm.
	 * (If {@link EuclideanUtils#NORM_EUCLIDEAN} was selected, then the supplied
	 * r should be the required Euclidean norm <b>squared</b>, since we switch it
	 * to {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED} internally).
	 */
	public int countPointsWithinR(int sampleIndex, double r, boolean allowEqualToR) {
		int count = 0;
		// Find where this node sits in the sorted array:
		int indexInSortedArray = indicesInSortedArray[sampleIndex];
		// Check the points with smaller data values first:
		for (int i = indexInSortedArray - 1; i >= 0; i--) {
			double theNorm = norm(originalDataSet[sampleIndex],
					originalDataSet[sortedArrayIndices[i]], normTypeToUse);
			if ((allowEqualToR  && (theNorm <= r)) ||
				(!allowEqualToR && (theNorm < r))) {
				count++;
				continue;
			}
			// Else no point checking further points
			break;
		}
		// Next check the points with larger data values:
		for (int i = indexInSortedArray + 1; i < numObservations; i++) {
			double theNorm = norm(originalDataSet[sampleIndex],
					originalDataSet[sortedArrayIndices[i]], normTypeToUse);
			if ((allowEqualToR  && (theNorm <= r)) ||
				(!allowEqualToR && (theNorm < r))) {
				count++;
				continue;
			}
			// Else no point checking further points
			break;
		}
		return count;
	}
	
	/**
	 * Count the number of points within norm r for a given
	 *  sample index in the data set. Nodes within dynCorrExclTime
	 *  of sampleIndex are excluded from the search.
	 * Nearest neighbour function to compare to r is the specified norm.
	 * (If {@link EuclideanUtils#NORM_EUCLIDEAN} was selected, then the supplied
	 * r should be the required Euclidean norm <b>squared</b>, since we switch it
	 * to {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED} internally).
	 */
	public int countPointsWithinR(int sampleIndex, double r,
			int dynCorrExclTime, boolean allowEqualToR) {
		int count = 0;
		// Find where this node sits in the sorted array:
		int indexInSortedArray = indicesInSortedArray[sampleIndex];
		// Check the points with smaller data values first:
		for (int i = indexInSortedArray - 1; i >= 0; i--) {
			if (Math.abs(sampleIndex - sortedArrayIndices[i]) <= dynCorrExclTime) {
				// Can't count this point, but keep checking:
				continue;
			}
			double theNorm = norm(originalDataSet[sampleIndex],
					originalDataSet[sortedArrayIndices[i]], normTypeToUse);
			if ((allowEqualToR  && (theNorm <= r)) ||
				(!allowEqualToR && (theNorm < r))) {
				count++;
				continue;
			}
			// Else no point checking further points
			break;
		}
		// Next check the points with larger data values:
		for (int i = indexInSortedArray + 1; i < numObservations; i++) {
			if (Math.abs(sampleIndex - sortedArrayIndices[i]) <= dynCorrExclTime) {
				// Can't count this point, but keep checking:
				continue;
			}
			double theNorm = norm(originalDataSet[sampleIndex],
					originalDataSet[sortedArrayIndices[i]], normTypeToUse);
			if ((allowEqualToR  && (theNorm <= r)) ||
				(!allowEqualToR && (theNorm < r))) {
				count++;
				continue;
			}
			// Else no point checking further points
			break;
		}
		return count;
	}
	
	/**
	 * Find the collection of points within norm r for a given
	 *  sample index in the data set. The node itself is 
	 *  excluded from the search.
	 * Nearest neighbour function to compare to r is the specified norm.
	 * (If {@link EuclideanUtils#NORM_EUCLIDEAN} was selected, then the supplied
	 * r should be the required Euclidean norm <b>squared</b>, since we switch it
	 * to {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED} internally).
	 */
	public Collection<NeighbourNodeData> findPointsWithinR(int sampleIndex, double r, boolean allowEqualToR) {
		Vector<NeighbourNodeData> pointsWithinR = new Vector<NeighbourNodeData>();
		// Find where this node sits in the sorted array:
		int indexInSortedArray = indicesInSortedArray[sampleIndex];
		// Check the points with smaller data values first:
		for (int i = indexInSortedArray - 1; i >= 0; i--) {
			double theNorm = norm(originalDataSet[sampleIndex],
					originalDataSet[sortedArrayIndices[i]], normTypeToUse);
			if ((allowEqualToR  && (theNorm <= r)) ||
				(!allowEqualToR && (theNorm < r))) {
				pointsWithinR.add(
					new NeighbourNodeData(sortedArrayIndices[i],
							new double[] {theNorm}, theNorm));
				continue;
			}
			// Else no point checking further points
			break;
		}
		// Next check the points with larger data values:
		for (int i = indexInSortedArray + 1; i < numObservations; i++) {
			double theNorm = norm(originalDataSet[sampleIndex],
					originalDataSet[sortedArrayIndices[i]], normTypeToUse);
			if ((allowEqualToR  && (theNorm <= r)) ||
				(!allowEqualToR && (theNorm < r))) {
				pointsWithinR.add(
						new NeighbourNodeData(sortedArrayIndices[i],
								new double[] {theNorm}, theNorm));
				continue;
			}
			// Else no point checking further points
			break;
		}
		return pointsWithinR;
	}
	
	/**
	 * Record the collection of points within norm r for a given
	 *  sample index in the data set. The node itself is 
	 *  excluded from the search.
	 * Nearest neighbour function to compare to r is the specified norm.
	 * (If {@link EuclideanUtils#NORM_EUCLIDEAN} was selected, then the supplied
	 * r should be the required Euclidean norm <b>squared</b>, since we switch it
	 * to {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED} internally).
	 * 
	 * <p>The recording of nearest neighbours is made within the isWithinR
	 *  and indicesWithinR arrays, which must be constructed before
	 *  calling this method, with length at or exceeding the total
	 *  number of data points. indicesWithinR is 
	 * </p> 
	 * 
	 */
	public void findPointsWithinR(
			int sampleIndex, double r, boolean allowEqualToR,
			boolean[] isWithinR, int[] indicesWithinR) {
		int indexInIndicesWithinR = 0;
		// Find where this node sits in the sorted array:
		int indexInSortedArray = indicesInSortedArray[sampleIndex];
		// Check the points with smaller data values first:
		for (int i = indexInSortedArray - 1; i >= 0; i--) {
			double theNorm = norm(originalDataSet[sampleIndex],
					originalDataSet[sortedArrayIndices[i]], normTypeToUse);
			if ((allowEqualToR  && (theNorm <= r)) ||
				(!allowEqualToR && (theNorm < r))) {
				isWithinR[sortedArrayIndices[i]] = true;
				indicesWithinR[indexInIndicesWithinR++] = sortedArrayIndices[i];
				continue;
			}
			// Else no point checking further points
			break;
		}
		// Next check the points with larger data values:
		for (int i = indexInSortedArray + 1; i < numObservations; i++) {
			double theNorm = norm(originalDataSet[sampleIndex],
					originalDataSet[sortedArrayIndices[i]], normTypeToUse);
			if ((allowEqualToR  && (theNorm <= r)) ||
				(!allowEqualToR && (theNorm < r))) {
				isWithinR[sortedArrayIndices[i]] = true;
				indicesWithinR[indexInIndicesWithinR++] = sortedArrayIndices[i];
				continue;
			}
			// Else no point checking further points
			break;
		}
		// Write the terminating integer into the indicesWithinR array:
		indicesWithinR[indexInIndicesWithinR++] = -1;
	}

	/**
	 * Record the collection of points within norm r for a given
	 *  sample index in the data set. The node itself is 
	 *  excluded from the search.
	 * Nearest neighbour function to compare to r is the specified norm.
	 * (If {@link EuclideanUtils#NORM_EUCLIDEAN} was selected, then the supplied
	 * r should be the required Euclidean norm <b>squared</b>, since we switch it
	 * to {@link EuclideanUtils#NORM_EUCLIDEAN_SQUARED} internally).
	 * 
	 * <p>The recording of nearest neighbours is made within the isWithinR
	 *  and indicesWithinR arrays, which must be constructed before
	 *  calling this method, with length at or exceeding the total
	 *  number of data points. indicesWithinR is 
	 * </p> 
	 * 
	 */
	public void findPointsWithinR(
			int sampleIndex, double r, int dynCorrExclTime,
			boolean allowEqualToR,
			boolean[] isWithinR, int[] indicesWithinR) {
		int indexInIndicesWithinR = 0;
		// Find where this node sits in the sorted array:
		int indexInSortedArray = indicesInSortedArray[sampleIndex];
		// Check the points with smaller data values first:
		for (int i = indexInSortedArray - 1; i >= 0; i--) {
			if (Math.abs(sampleIndex - sortedArrayIndices[i]) <= dynCorrExclTime) {
				// Can't count this point, but keep checking:
				continue;
			}
			double theNorm = norm(originalDataSet[sampleIndex],
					originalDataSet[sortedArrayIndices[i]], normTypeToUse);
			if ((allowEqualToR  && (theNorm <= r)) ||
				(!allowEqualToR && (theNorm < r))) {
				isWithinR[sortedArrayIndices[i]] = true;
				indicesWithinR[indexInIndicesWithinR++] = sortedArrayIndices[i];
				continue;
			}
			// Else no point checking further points
			break;
		}
		// Next check the points with larger data values:
		for (int i = indexInSortedArray + 1; i < numObservations; i++) {
			if (Math.abs(sampleIndex - sortedArrayIndices[i]) <= dynCorrExclTime) {
				// Can't count this point, but keep checking:
				continue;
			}
			double theNorm = norm(originalDataSet[sampleIndex],
					originalDataSet[sortedArrayIndices[i]], normTypeToUse);
			if ((allowEqualToR  && (theNorm <= r)) ||
				(!allowEqualToR && (theNorm < r))) {
				isWithinR[sortedArrayIndices[i]] = true;
				indicesWithinR[indexInIndicesWithinR++] = sortedArrayIndices[i];
				continue;
			}
			// Else no point checking further points
			break;
		}
		// Write the terminating integer into the indicesWithinR array:
		indicesWithinR[indexInIndicesWithinR++] = -1;
	}
	
	@Override
	public int countPointsWithinR(int sampleIndex, double r, boolean allowEqualToR,
			boolean[] additionalCriteria) {
		
		int count = 0;
		// Find where this node sits in the sorted array:
		int indexInSortedArray = indicesInSortedArray[sampleIndex];
		// Check the points with smaller data values first:
		for (int i = indexInSortedArray - 1; i >= 0; i--) {
			if (!additionalCriteria[sortedArrayIndices[i]]) {
				// This point failed the additional criteria,
				//  so skip checking it.
				// We check this first, even though it's the norm 
				//  that really determines whether we need to continue
				//  checking or not. In some circumstances this may be
				//  slower, but for our main application - KSG conditional MI -
				//  this should be faster.
				continue;
			}
			double theNorm = norm(originalDataSet[sampleIndex],
					originalDataSet[sortedArrayIndices[i]], normTypeToUse);
			if ((allowEqualToR  && (theNorm <= r)) ||
				(!allowEqualToR && (theNorm < r))) {
				count++;
				continue;
			}
			// Else no point checking further points
			break;
		}
		// Next check the points with larger data values:
		for (int i = indexInSortedArray + 1; i < numObservations; i++) {
			if (!additionalCriteria[sortedArrayIndices[i]]) {
				// This point failed the additional criteria,
				//  so skip checking it.
				// We check this first, even though it's the norm 
				//  that really determines whether we need to continue
				//  checking or not. In some circumstances this may be
				//  slower, but for our main application - KSG conditional MI -
				//  this should be faster.
				continue;
			}
			double theNorm = norm(originalDataSet[sampleIndex],
					originalDataSet[sortedArrayIndices[i]], normTypeToUse);
			if ((allowEqualToR  && (theNorm <= r)) ||
				(!allowEqualToR && (theNorm < r))) {
				count++;
				continue;
			}
			// Else no point checking further points
			break;
		}
		return count;
	}
}
