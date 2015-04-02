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

/**
 * 
 * A set of utilities for manipulating vectors in Euclidean space.
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class EuclideanUtils {

	public static int MAX_TIMESTEPS_FOR_FAST_DISTANCE = 2000;
	
	public static final int NORM_EUCLIDEAN = 0;
	public static final String NORM_EUCLIDEAN_STRING = "EUCLIDEAN";
	public static final int NORM_EUCLIDEAN_NORMALISED = 1;
	public static final String NORM_EUCLIDEAN_NORMALISED_STRING = "EUCLIDEAN_NORMALISED";
	public static final int NORM_MAX_NORM = 2;
	public static final String NORM_MAX_NORM_STRING = "MAX_NORM";
	public static final int NORM_EUCLIDEAN_SQUARED = 3;
	public static final String NORM_EUCLIDEAN_SQUARED_STRING = "EUCLIDEAN_SQUARED";
	// Track which norm we should use here
	private int normToUse = 0;

	/**
	 * Construct a EuclideanUtils object, to take norms of the given type
	 * 
	 * @param normToUse norm type, one of
	 * {@link #NORM_EUCLIDEAN},
	 * {@link #NORM_EUCLIDEAN_NORMALISED},
	 * {@link #NORM_MAX_NORM}
	 *  or {@link #NORM_MAX_NORM}
	 */
	public EuclideanUtils(int normToUse) {
		setNormToUse(normToUse);
	}

	/**
	 * Construct a EuclideanUtils object, to take norms of the given type
	 * 
	 * @param normToUse norm type, one of
	 * {@link #NORM_EUCLIDEAN_STRING},
	 * {@link #NORM_EUCLIDEAN_NORMALISED_STRING},
	 * {@link #NORM_MAX_NORM_STRING}
	 *  or {@link #NORM_MAX_NORM_STRING}
	 */
	public EuclideanUtils(String normToUse) {
		setNormToUse(normToUse);
	}

	public static double[] computeMinEuclideanDistances(double[][] observations) {
		if (observations.length <= MAX_TIMESTEPS_FOR_FAST_DISTANCE) {
			return EuclideanUtils.computeMinEuclideanDistancesFast(observations);
		} else {
			return computeMinEuclideanDistancesNaive(observations);
		}
	}
	
	
	/**
	 * Naive method for computing minimum distance - slower but needs less memory.
	 * Made public for debugging only. O(d.n^2) speed
	 * 
	 * @param observations
	 * @return
	 */
	public static double[] computeMinEuclideanDistancesNaive(double[][] observations) {
		int numObservations = observations.length;
		int dimensions = observations[0].length;
		double[] distances = new double[numObservations];
		for (int t = 0; t < numObservations; t++) {
			double minDistance = Double.POSITIVE_INFINITY;
			for (int t2 = 0; t2 < numObservations; t2++) {
				if (t == t2) {
					continue;
				}
				double thisDistance = 0.0;
				for (int d = 0; (d < dimensions) && (thisDistance < minDistance); d++) {
					double distanceOnThisVar = (observations[t][d] - observations[t2][d]);
					thisDistance += distanceOnThisVar * distanceOnThisVar;
				}
				// Now we need to sqrt the distance sum
				thisDistance = Math.sqrt(thisDistance);
				// Now check if this is a lower distance
				if (thisDistance < minDistance) {
					minDistance = thisDistance;
				}
			}
			distances[t] = minDistance;
		}
		return distances;
	}

	/**
	 * Return the minimum Euclidean distance from each point to any other observation.
	 * Computes this faster than using naive computation.
	 * 
	 * Exposed as a public method for debugging purposes only.
	 * 
	 * @param observations
	 * @return
	 */
	public static double[] computeMinEuclideanDistancesFast(double[][] observations) {
		
		int dimensions  = observations[0].length;
		
		int timeSteps = observations.length;
		// Hold the sqr distance from index1 to index2 ...
		double[][] sqrDistance = new double[timeSteps][timeSteps];
		// ... computed over this many of the variables so far
		int[][] addedInUpToVariable = new int[timeSteps][timeSteps];
		double[] minDistance = new double[timeSteps];
		
		for (int t1 = 0; t1 < timeSteps; t1++) {
			// Current minimum distance from this index to another point:
			double minSqrDistance = Double.POSITIVE_INFINITY;

			// First grab the minimum distance from nodes for which the distance might
			//  have already been measured
			for (int t2 = 0; t2 < t1; t2++) {
				if (addedInUpToVariable[t2][t1] == dimensions) {
					// We have previously computed this distance from t2 to t1
					sqrDistance[t1][t2] = sqrDistance[t2][t1];
					// unnecessary, since we won't be looking at [t1][t2] later:
					addedInUpToVariable[t1][t2] = dimensions;
					if (sqrDistance[t1][t2] < minSqrDistance) {
						minSqrDistance = sqrDistance[t1][t2];
					}
				}
			}
			// Now check the previously considered source nodes which didn't have their full distance
			//  computed in case we need to compute them
			for (int t2 = 0; t2 < t1; t2++) {
				if (addedInUpToVariable[t2][t1] != dimensions) {
					// We have not finished computing this distance from t1
					addedInUpToVariable[t1][t2] = addedInUpToVariable[t2][t1];
					sqrDistance[t1][t2] = sqrDistance[t2][t1];
					for (; (sqrDistance[t1][t2] < minSqrDistance) &&
							(addedInUpToVariable[t1][t2] < dimensions);
							addedInUpToVariable[t1][t2]++) {
						double distOnThisVar = observations[t1][addedInUpToVariable[t1][t2]] -
												observations[t2][addedInUpToVariable[t1][t2]]; 
						sqrDistance[t1][t2] += distOnThisVar * distOnThisVar;
					}
					if (sqrDistance[t1][t2] < minSqrDistance) {
						// we finished the calculation and t2 is now the closest observation to t1
						minSqrDistance = sqrDistance[t1][t2];
					}
				}
			}
			// Now check any source nodes t2 for which there is no chance we've looked at the
			//  the distance back to t1 yet
			for (int t2 = t1 + 1; t2 < timeSteps; t2++) {
				for (; (sqrDistance[t1][t2] < minSqrDistance) &&
						(addedInUpToVariable[t1][t2] < dimensions);
						addedInUpToVariable[t1][t2]++) {
					double distOnThisVar = observations[t1][addedInUpToVariable[t1][t2]] -
											observations[t2][addedInUpToVariable[t1][t2]]; 
					sqrDistance[t1][t2] += distOnThisVar * distOnThisVar;
				}
				if (sqrDistance[t1][t2] < minSqrDistance) {
					// we finished the calculation and  t2 is now the closest observation to t1
					minSqrDistance = sqrDistance[t1][t2];
				}
			}
			
			minDistance[t1] = Math.sqrt(minSqrDistance);
		}
		return minDistance;
	}

	/**
	 * Return the max norm out of the two norms (x1:x2) and (y1:y2),
	 *  using the configured norm type
	 * 
	 * @param x1
	 * @param y1
	 * @param x2
	 * @param y2
	 * @return the max of the two norms
	 */
	public double maxJointSpaceNorm(double[] x1, double[] y1,
			   double[] x2, double[] y2) {
		return Math.max(norm(x1, x2), norm(y1,y2));
	}

	/**
	 * Computing the configured norm between vectors x1 and x2.
	 * 
	 * @param x1
	 * @param x2
	 * @return
	 */
	public double norm(double[] x1, double[] x2) {
		switch (normToUse) {
		case NORM_EUCLIDEAN_NORMALISED:
			return euclideanNorm(x1, x2) / Math.sqrt(x1.length);
		case NORM_MAX_NORM:
			return maxNorm(x1, x2);
		case NORM_EUCLIDEAN_SQUARED:
			return euclideanNormSquared(x1, x2);
		case NORM_EUCLIDEAN:
		default:
			return euclideanNorm(x1, x2);
		}
	}
	
	/**
	 * Computing the configured norm between vectors x1 and x2; if 
	 *  it becomes clear that norm will be larger than limit,
	 *  then return Double.POSITIVE_INFINITY immediately.
	 * 
	 * @param x1
	 * @param x2
	 * @param limit
	 * @return
	 */
	public double normWithAbort(double[] x1, double[] x2, double limit) {
		switch (normToUse) {
		case NORM_EUCLIDEAN_NORMALISED:
			return euclideanNormWithAbort(x1, x2, limit) / Math.sqrt(x1.length);
		case NORM_MAX_NORM:
			return maxNormWithAbort(x1, x2, limit);
		case NORM_EUCLIDEAN_SQUARED:
			return euclideanNormSquaredWithAbort(x1, x2, limit);
		case NORM_EUCLIDEAN:
		default:
			return euclideanNormWithAbort(x1, x2, limit);
		}
	}
	
	/**
	 * Computing the norm as the Euclidean norm.
	 * 
	 * @param x1
	 * @param x2
	 * @return
	 */
	public static double euclideanNorm(double[] x1, double[] x2) {
		double distance = 0.0;
		for (int d = 0; d < x1.length; d++) {
			double difference = x1[d] - x2[d];
			distance += difference * difference;
		}
		return Math.sqrt(distance);
	}
	
	/**
	 * Computing the norm as the Euclidean norm; if 
	 *  it becomes clear that norm will be larger than limit,
	 *  then return Double.POSITIVE_INFINITY immediately.
	 * 
	 * @param x1
	 * @param x2
	 * @param limit
	 * @return
	 */
	public static double euclideanNormWithAbort(double[] x1, double[] x2, double limit) {
		double distance = 0.0;
		limit *= limit;
		for (int d = 0; d < x1.length; d++) {
			double difference = x1[d] - x2[d];
			distance += difference * difference;
			if (distance > limit) {
				return Double.POSITIVE_INFINITY;
			}
		}
		return Math.sqrt(distance);
	}

	/**
	 * Computing the norm as the Euclidean norm squared
	 *  (i.e. avoids taking the square root).
	 * 
	 * @param x1
	 * @param x2
	 * @return
	 */
	public static double euclideanNormSquared(double[] x1, double[] x2) {
		double distance = 0.0;
		for (int d = 0; d < x1.length; d++) {
			double difference = x1[d] - x2[d];
			distance += difference * difference;
		}
		return distance;
	}
	
	/**
	 * Computing the norm as the Euclidean norm squared
	 *  (i.e. avoids taking the square root); if 
	 *  it becomes clear that norm will be larger than limit,
	 *  then return Double.POSITIVE_INFINITY immediately.
	 * 
	 * @param x1 vector 1
	 * @param x2 vector 2
	 * @return
	 */
	public static double euclideanNormSquaredWithAbort(
			double[] x1, double[] x2, double limit) {
		double distance = 0.0;
		for (int d = 0; d < x1.length; d++) {
			double difference = x1[d] - x2[d];
			distance += difference * difference;
			if (distance > limit) {
					return Double.POSITIVE_INFINITY;
			}
		}
		return distance;
	}
	
	/**
	 * Computing the norm as the Max norm.
	 * 
	 * @param x1
	 * @param x2
	 * @return
	 */
	public static double maxNorm(double[] x1, double[] x2) {
		double distance = 0.0;
		for (int d = 0; d < x1.length; d++) {
			double difference = x1[d] - x2[d];
			// Take the abs
			if (difference < 0) {
				difference = -difference;
			}
			if (difference > distance) {
				distance = difference;
			}
		}
		return distance;
	}

	/**
	 * Computing the norm as the Max norm; if 
	 *  it becomes clear that norm will be larger than limit,
	 *  then return Double.POSITIVE_INFINITY immediately.
	 * 
	 * @param x1 vector 1
	 * @param x2 vector 2
	 * @param limit
	 * @return
	 */
	public static double maxNormWithAbort(double[] x1, double[] x2, double limit) {
		double distance = 0.0;
		for (int d = 0; d < x1.length; d++) {
			double difference = x1[d] - x2[d];
			// Take the abs
			if (difference < 0) {
				difference = -difference;
			}
			if (difference > distance) {
				if (difference > limit) {
					return Double.POSITIVE_INFINITY;
				}
				distance = difference;
			}
		}
		return distance;
	}
	
	/**
	 * Compute the x and y configured norms of all other points from
	 *  the data points at time step t.
	 * Puts norms of t from itself as infinity, which is useful
	 *  when counting the number of points closer than epsilon say.
	 * 
	 * @param mvTimeSeries1
	 * @param mvTimeSeries2
	 * @return
	 */
	public double[][] computeNorms(double[][] mvTimeSeries1,
			double[][] mvTimeSeries2, int t) {
		
		int timeSteps = mvTimeSeries1.length;
		double[][] norms = new double[timeSteps][2];
		for (int t2 = 0; t2 < timeSteps; t2++) {
			if (t2 == t) {
				norms[t2][0] = Double.POSITIVE_INFINITY;
				norms[t2][1] = Double.POSITIVE_INFINITY;
				continue;
			}
			// Compute norm in first direction
			norms[t2][0] = norm(mvTimeSeries1[t], mvTimeSeries1[t2]);
			// Compute norm in second direction
			norms[t2][1] = norm(mvTimeSeries2[t], mvTimeSeries2[t2]);
		}
		return norms;
	}
	
	/**
	 * Compute the x, y and z norms of all other points from
	 *  the data points at time step t.
	 * Puts norms of t from itself as infinity, which is useful
	 *  when counting the number of points closer than epsilon say.
	 * 
	 * @param mvTimeSeries1
	 * @param mvTimeSeries2
	 * @param mvTimeSeries3
	 * @return
	 */
	public double[][] computeNorms(double[][] mvTimeSeries1,
			double[][] mvTimeSeries2, double[][] mvTimeSeries3, int t) {
		
		int timeSteps = mvTimeSeries1.length;
		double[][] norms = new double[timeSteps][3];
		for (int t2 = 0; t2 < timeSteps; t2++) {
			if (t2 == t) {
				norms[t2][0] = Double.POSITIVE_INFINITY;
				norms[t2][1] = Double.POSITIVE_INFINITY;
				norms[t2][2] = Double.POSITIVE_INFINITY;
				continue;
			}
			// Compute norm in first direction
			norms[t2][0] = norm(mvTimeSeries1[t], mvTimeSeries1[t2]);
			// Compute norm in second direction
			norms[t2][1] = norm(mvTimeSeries2[t], mvTimeSeries2[t2]);
			// Compute norm in third direction
			norms[t2][2] = norm(mvTimeSeries3[t], mvTimeSeries3[t2]);
		}
		return norms;
	}

	/**
	 * Compute the norms for each marginal variable for all other points from
	 *  the data points at time step t.
	 * Puts norms of t from itself as infinity, which is useful
	 *  when counting the number of points closer than epsilon say.
	 * 
	 * @param mvTimeSeries
	 * @return
	 */
	public static double[][] computeNorms(double[][] mvTimeSeries, int t) {
		
		int timeSteps = mvTimeSeries.length;
		int variables = mvTimeSeries[0].length;
		
		double[][] norms = new double[timeSteps][variables];
		for (int t2 = 0; t2 < timeSteps; t2++) {
			if (t2 == t) {
				for (int v = 0; v < variables; v++) {
					norms[t2][v] = Double.POSITIVE_INFINITY;
				}
				continue;
			}
			for (int v = 0; v < variables; v++) {
				norms[t2][v] = Math.abs(mvTimeSeries[t][v] - mvTimeSeries[t2][v]);
			}
		}
		return norms;
	}

	/**
	 * Sets which type of norm will be used by calls to norm()
	 * 
	 * @param normType
	 */
	public void setNormToUse(int normType) {
		this.normToUse = normType;
	}

	/**
	 * Sets which type of norm will be used by calls to norm()
	 * 
	 * @param normType
	 */
	public void setNormToUse(String normType) {
		if (normType.equalsIgnoreCase(NORM_EUCLIDEAN_NORMALISED_STRING)) {
			normToUse = NORM_EUCLIDEAN_NORMALISED;
		} else if (normType.equalsIgnoreCase(NORM_MAX_NORM_STRING)) {
			normToUse = NORM_MAX_NORM;
		} else if (normType.equalsIgnoreCase(NORM_EUCLIDEAN_SQUARED_STRING)) {
			normToUse = NORM_EUCLIDEAN_SQUARED;
		} else {
			normToUse = NORM_EUCLIDEAN;
		}
	}
	
	public int getNormInUse() {
		return normToUse;
	}

	public String getNormInUseString() {
		switch (normToUse) {
		case NORM_EUCLIDEAN_NORMALISED:
			return NORM_EUCLIDEAN_NORMALISED_STRING;
		case NORM_MAX_NORM:
			return NORM_MAX_NORM_STRING;
		case NORM_EUCLIDEAN_SQUARED:
			return NORM_EUCLIDEAN_SQUARED_STRING;
		default:
		case NORM_EUCLIDEAN:
			return NORM_EUCLIDEAN_STRING;
		}
	}
}
