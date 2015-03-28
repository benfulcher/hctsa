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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Random;
import java.util.Vector;
import java.util.Hashtable;
import java.util.Arrays;

/**
 * Utility to generate arrays of random variables
 * 
 * TODO I think I may need to revisit whether I reuse the objects added to the hashtable;
 * I don't think we should be doing this.
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class RandomGenerator {

	Random random;
	
	public RandomGenerator() {
		random = new Random();
	}

	public void setSeed(long seed) {
		random.setSeed(seed);
	}
	
	public double[] generateRandomData(int length){
		double[] data = new double[length];
		for (int i = 0; i < length; i++) {
			data[i] = random.nextDouble();
		}
		return data;
	}

	/**
	 * Generate an array of random data on the interval [min .. max)
	 * 
	 * @param length
	 * @param min
	 * @param max
	 * @return
	 */
	public double[] generateRandomData(int length, double min, double max){
		double[] data = new double[length];
		for (int i = 0; i < length; i++) {
			data[i] = min + random.nextDouble() * (max - min);
		}
		return data;
	}

	public double[][] generateRandomData(int length, int dimenions){
		double[][] data = new double[length][dimenions];
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < dimenions; j++) {
				data[i][j] = random.nextDouble();
			}
		}
		return data;
	}
	
	/**
	 * <p>
	 * Generate <i>length</i> random ints, between the values 0..(<i>cap</i>-1)
	 * </p>
	 * 
	 * @param length
	 * @param cap
	 * @return
	 */
	public int[] generateDistinctRandomInts(int length, int cap){
		int[] data = new int[length];
		boolean[] used = new boolean[cap];
		for (int i = 0; i < length; i++) {
			int nextAttempt;
			// Select an int we haven't used yet:
			for (nextAttempt = random.nextInt(cap); used[nextAttempt]; nextAttempt = random.nextInt(cap)) {
				// Select the next attempt
			}
			data[i] = nextAttempt;
			used[nextAttempt] = true;
		}
		return data;
	}

	public double[] generateNormalData(int length, double mean, double std){
		double[] data = new double[length];
		for (int i = 0; i < length; i++) {
			data[i] = random.nextGaussian()*std + mean;
		}
		return data;
	}

	public double[][] generateNormalData(int length, int dimensions,
			double mean, double std){
		double[][] data = new double[length][dimensions];
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < dimensions; j++) {
				data[i][j] = random.nextGaussian()*std + mean;
			}
		}
		return data;
	}

	/**
	 * <p>Generate bivariate Gaussian series with the given covariance.
	 * See http://mathworld.wolfram.com/BivariateNormalDistribution.html</p>
	 * <p>
	 * If we have two normal distributions x1 and x2, we can define</br>
	 *   <ul>
	 *   	<li>y1 = mean1 + sigma11*x1 + sigma12*x2</li>
	 *   	<li>y2 = mean2 + sigma21*x1 + sigma22*x2</li>
	 *   </ul>
	 * which are Gaussian distributed with:
	 * 	<ul>
	 *  	<li>means (mean1,mean2),</li>
	 *  	<li>variances (sigma11^2+sigma12^2, sigma21^2+sigma22^2), and</li>
	 *  	<li>covariance sigma11*sigma21 + sigma12*sigma22</li>
	 *  </ul>
	 * So to generate a bivariate series with a desired covariance, means and stds,
	 * we set sigma12=0, giving sigma11 = std1, solve for sigma21 from the covariance,
	 * and solve for sigma22 from the std2.
	 * </p> 
	 * 
	 * @param length
	 * @param mean1
	 * @param std1
	 * @param mean2
	 * @param std2
	 * @param covariance
	 * @return a time series with two variables: first index is time step, second index is variable number
	 */
	public double[][] generateBivariateNormalData(int length,
			double mean1, double std1,
			double mean2, double std2,
			double covariance){

		double[][] data = new double[length][2];
		double sigma21 = covariance / std1;
		double sigma22 = Math.sqrt(std2*std2 - sigma21*sigma21);
		for (int i = 0; i < length; i++) {
			double x1 = random.nextGaussian();
			double x2 = random.nextGaussian();
			data[i][0] = mean1 + std1*x1;
			data[i][1] = mean2 + sigma21*x1 + sigma22*x2;
		}
		return data;
	}
	
	/**
	 * Generate a set of covariant gaussians, with the given means and covariances.
	 * 
	 * @param length Number of time steps (samples) generated
	 * @param dimensions Number of gaussians
	 * @param means Means of the generated gaussians
	 * @param componentDependencies Underlying depencence matrix A between the generated gaussians.
	 * 		Covariance C = A * A^T
	 * @return
	 */
	public double[][] generateCovariantGaussians(int length, int dimensions,
			double[] means, double[][] componentDependencies) {
		double[][] data = new double[length][dimensions];
		for (int t = 0; t < length; t++) {
			// Generate the underlying random values for this time step
			double[] x = generateNormalData(dimensions, 0, 1);
			for (int d = 0; d < dimensions; d++) {
				data[t][d] = means[d];
				// Combine the underlying random values for variable d
				for (int d2 = 0; d2 < dimensions; d2++) {
					data[t][d] += componentDependencies[d][d2] * x[d2];
				}
			}
		}
		return data;
	}
	
	/**
	 * Generate an array of random integers in the range 0 .. cap-1.
	 * 
	 * @param length length of array to return
	 * @param cap number of distinct values to choose from (0..cap-1)
	 * @return array of random integers
	 */
	public int[] generateRandomInts(int length, int cap) {
		int[] data = new int[length];
		for (int i = 0; i < length; i++) {
			data[i] = random.nextInt(cap);
		}
		return data;
	}
	
	/**
	 * Generate a multidimensional array of random integers in the range 0 .. cap-1.
	 * 
	 * @param rows rows of array to return
	 * @param columns columns of array to return
	 * @param cap number of distinct values to choose from (0..cap-1)
	 * @return array of random integers
	 */
	public int[][] generateRandomInts(int rows, int columns, int cap) {
		int[][] data = new int[rows][columns];
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				data[r][c] = random.nextInt(cap);
			}
		}
		return data;
	}
	
	/**
	 * Generate (up to) numSets distinct sets of p distinct values in [0..n-1] 
	 * Done using random guesses as this is designed for high dimension n
	 *  where its highly unlikely we repeat a set (though this is checked)
	 * 
	 * @param n
	 * @param p
	 * @param maxNumSets
	 * @return
	 */
	public int[][] generateDistinctRandomSets(int n, int p, int maxNumSets) {
		
		// Check what is the max possible number of sets we
		//  could generate:
		int maxPossibleNumSets = 0;
		try {
			maxPossibleNumSets = MathsUtils.numOfSets(n, p);
			if (maxNumSets > maxPossibleNumSets) {
				// Best limit maxNumSets
				maxNumSets = maxPossibleNumSets;
				// We want to generate all possible sets
				return generateAllDistinctSets(n, p);
			}
		} catch (Exception e) {
			// n choose p blew Integer.MAX_INT
			// therefore there is no way maxNumSets is larger than it
		}
		
		int[][] sets = new int[maxNumSets][p];
		// Pool of available choices:
		Vector<Integer> availableChoices = new Vector<Integer>();
		for (int i = 0; i < n; i++) {
			availableChoices.add(new Integer(i));
		}
		// Pool of choices already used this turn
		Vector<Integer> thisSet = new Vector<Integer>();
		// Hashtable tracking sets we've already chosen
		Hashtable<Vector<Integer>,Integer> chosenSets =
			new Hashtable<Vector<Integer>,Integer>();
		for (int s = 0; s < maxNumSets; s++) {
			// Select set s:
			for (;;) {
				// Try to get a new unique set
				// Reset the pool of choices ready to choose this set
				availableChoices.addAll(thisSet);
				thisSet.clear();
				// System.out.println("Available: " + availableChoices);
				for (int q = 0; q < p; q++) {
					// Select the qth index from the available pool to use here:
					int randIndex = random.nextInt(n - q);
					// Find out what number this corresponds to, and write it in:
					Integer nextSelection = availableChoices.remove(randIndex);
					sets[s][q] = nextSelection.intValue();
				}
				// Track the chosen integers in order to avoid duplicates
				Arrays.sort(sets[s]);
				for (int q = 0; q < p; q++) {
					// And track it in thisSet for hashing and 
					//  adding back to the pool
					thisSet.add(new Integer(sets[s][q]));
				}
				if (chosenSets.get(thisSet) == null) {
					// We haven't included this set yet:
					chosenSets.put(thisSet, new Integer(0));
					// System.out.println(" Chosen: " + thisSet);
					break;
				}
				// else we have already added this set,
				//  so we need to try again.
				// System.out.println(" Attempted: " + thisSet);
				// System.out.println(" Need to try again");
			}
		}
		return sets;
	}
	
	/**
	 * Generate exactly N random sets of p numbers from [0 .. n-1],
	 *  allowing repeats if nCp < N.
	 * 
	 * @param n
	 * @param p
	 * @param N
	 * @return
	 */
	public int[][] generateNRandomSets(int n, int p,
			int N) {
		int[][] distinctSets = generateDistinctRandomSets(n, p, N);
		if (distinctSets.length == N) {
			// Fine - return it
			return distinctSets;
		} else if (distinctSets.length > N) {
			// Error condition - generateDistinctRandomSets
			//  should not do this
			throw new RuntimeException(
				"generateDistinctRandomSets generated more than " +
				N + " distinct sets when asked for " + N + 
				"; note n=" + n + " p=" + p);
		}
		// Else we now scale up the available distinct rows
		//  to fill the whole N required sets
		int[][] randomSets = new int[N][];
		for (int i = 0; i < N; i++) {
			// Select one of the distinct sets at random:
			randomSets[i] = distinctSets[random.nextInt(distinctSets.length)];
		}
		return randomSets;
	}

	/**
	 * Generate (up to) numSets distinct sets of p distinct values in [0..n-1].
	 * Done using random guesses as this is designed for high dimension n
	 *  where its highly unlikely we repeat a set (though this is checked).
	 * Here we avoid overlapping the sets with any elements of the corresponding
	 *  row of setsToAvoidOverlapWith
	 * 
	 * @param n
	 * @param p
	 * @param maxNumSets
	 * @param setsToAvoidOverlapWith
	 * @return
	 */
	private int[][] generateDistinctRandomSets(int n, int p, int maxNumSets, int[][] setsToAvoidOverlapWith) {
		
		// TODO Write this
		
		if (true) {
			throw new RuntimeException("Not implemented yet");
		}
		
		// Not sure whether we need to make this call beforehand or not - have a think about
		
		
		// Check what is the max possible number of sets we
		//  could generate:
		int maxPossibleNumSets = 0;
		try {
			maxPossibleNumSets = MathsUtils.numOfSets(n, p);
			if (maxNumSets > maxPossibleNumSets) {
				// Best limit maxNumSets
				maxNumSets = maxPossibleNumSets;
				// We want to generate all possible sets
				return generateAllDistinctSets(n, p);
			}
		} catch (Exception e) {
			// n choose p blew Integer.MAX_INT
			// therefore there is no way maxNumSets is larger than it
		}
		
		int[][] sets = new int[maxNumSets][p];
		// Pool of available choices:
		Vector<Integer> availableChoices = new Vector<Integer>();
		for (int i = 0; i < n; i++) {
			availableChoices.add(new Integer(i));
		}
		// Pool of choices already used this turn
		Vector<Integer> thisSet = new Vector<Integer>();
		// Hashtable tracking sets we've already chosen
		Hashtable<Vector<Integer>,Integer> chosenSets =
			new Hashtable<Vector<Integer>,Integer>();
		for (int s = 0; s < maxNumSets; s++) {
			// Select set s:
			for (;;) {
				// Try to get a new unique set
				// Reset the pool of choices ready to choose this set
				availableChoices.addAll(thisSet);
				// And remove the integers already chosen in the corresponding row of setsToAvoidOverlapWith:
				if (setsToAvoidOverlapWith != null) {
					for (int i = 0; i < p; i++) {
						availableChoices.remove(new Integer(setsToAvoidOverlapWith[s][p]));
					}
				}
				thisSet.clear();
				// System.out.println("Available: " + availableChoices);
				for (int q = 0; q < p; q++) {
					// Select the qth index from the available pool to use here:
					int randIndex = random.nextInt(n - q);
					// Find out what number this corresponds to, and write it in:
					Integer nextSelection = availableChoices.remove(randIndex);
					sets[s][q] = nextSelection.intValue();
				}
				// Track the chosen integers in order to avoid duplicates
				Arrays.sort(sets[s]);
				for (int q = 0; q < p; q++) {
					// And track it in thisSet for hashing and 
					//  adding back to the pool
					thisSet.add(new Integer(sets[s][q]));
				}
				if (chosenSets.get(thisSet) == null) {
					// We haven't included this set yet.
					chosenSets.put(thisSet, new Integer(0));
					// System.out.println(" Chosen: " + thisSet);
					break;
				}
				// else we have already added this set,
				//  so we need to try again.
				// System.out.println(" Attempted: " + thisSet);
				// System.out.println(" Need to try again");
			}
		}
		return sets;
	}

	/**
	 * Generate exactly N sets of p integers chosen from [0..n-1].
	 * Make sure that set i does not have any integers overlapping with set i
	 *  from setsToAvoidOverlapWith.
	 * 
	 * @param n
	 * @param p
	 * @param N
	 * @param setsToAvoidOverlapWith
	 * @return
	 */
	public int[][] generateNRandomSetsNoOverlap(int n, int p,
			int N, int[][] setsToAvoidOverlapWith) {
		// We generate N random sets of p numbers from n,
		//  allowing repeats if nCp < N.
		int[][] distinctSets = generateDistinctRandomSets(n, p, N, setsToAvoidOverlapWith);
		if (distinctSets.length == N) {
			// Fine - return it
			return distinctSets;
		} else if (distinctSets.length > N) {
			// Error condition - generateDistinctRandomSets
			//  should not do this
			throw new RuntimeException(
				"generateDistinctRandomSets generated more than " +
				N + " distinct sets when asked for " + N + 
				"; note n=" + n + " p=" + p);
		}
		// Else we now scale up the available distinct rows
		//  to fill the whole N required sets
		int[][] randomSets = new int[N][];
		for (int i = 0; i < N; i++) {
			// Select one of the distinct sets at random:
			randomSets[i] = distinctSets[random.nextInt(distinctSets.length)];
		}
		return randomSets;
	}

	/**
	 * Generate all nCp sets (assuming this doesn't blow our memory
	 * 
	 * @param n
	 * @param p
	 * @return
	 */
	public int[][] generateAllDistinctSets(int n, int p) throws Exception {
		int numSets;
		try {
			numSets = MathsUtils.numOfSets(n, p);
		} catch (Exception e) {
			// n choose p blew Integer.MAX_INT
			throw new Exception("nCp too large");
		}
		// allocate space for the distinct sets
		int[][] sets = new int[numSets][p];
		int[] workingSet = new int[p];
		addToDistinctSets(sets, n, p, 0, workingSet, 0, 0);
		return sets;
	}
	
	/**
	 * Using the workingSet which is filled up to (but not including)
	 *  fromIndex, add new distinct sets of nCp to the sets matrix,
	 *  from the setNumber index
	 * 
	 * @param sets
	 * @param n
	 * @param p
	 * @param setNumber
	 * @param workingSet
	 * @param fromIndex
	 */
	protected int addToDistinctSets(int[][] sets, int n, int p, int setNumber, 
			int[] workingSet, int fromIndex, int selectFrom) {
		if (fromIndex == p) {
			// The workingSet is ready to go, so copy it in
			// MatrixUtils.printArray(System.out, workingSet);
			System.arraycopy(workingSet, 0, sets[setNumber], 0, p);
			setNumber++;
		} else {
			// Add to the working set and pass it on:
			for (int c = selectFrom; c < n; c++) {
				workingSet[fromIndex] = c;
				setNumber = addToDistinctSets(sets, n, p, setNumber,
						workingSet, fromIndex + 1, c + 1);
			}
		}
		return setNumber;
	}
	
	/**
	 * Generate numberOfPerturbations perturbations of [0..n-1]
	 * 
	 * @param n
	 * @param numberOfPerturbations
	 * @return an array of dimensions [numberOfPerturbations][n], with each row
	 *  being one perturbation of the elements
	 */
	public int[][] generateDistinctRandomPerturbations(int n, int numberOfPerturbations) {
		// Check what is the max possible number of perturbations we
		//  could generate:
		int maxPossibleNumPerturbations = 0;
		try {
			maxPossibleNumPerturbations = MathsUtils.factorialCheckBounds(n);
			if (numberOfPerturbations > maxPossibleNumPerturbations) {
				// Best limit maxNumSets
				numberOfPerturbations = maxPossibleNumPerturbations;
				// We want to generate all possible sets
				return generateAllDistinctPerturbations(n);
			}
		} catch (Exception e) {
			// n! blew Integer.MAX_INT
			// therefore there is no way numberOfPerturbations is larger than it
		}
		
		int[][] sets = new int[numberOfPerturbations][n];
		// Pool of available choices:
		Vector<Integer> availableChoices = new Vector<Integer>();
		for (int i = 0; i < n; i++) {
			availableChoices.add(new Integer(i));
		}
		// Pool of choices already used this turn
		Vector<Integer> thisSet = new Vector<Integer>();
		// Hashtable tracking sets we've already chosen
		Hashtable<Vector<Integer>,Integer> chosenSets =
			new Hashtable<Vector<Integer>,Integer>();
		for (int s = 0; s < numberOfPerturbations; s++) {
			// Select set s:
			for (;;) {
				// Try to get a new unique set
				// Reset the pool of choices ready to choose this set
				availableChoices.addAll(thisSet);
				thisSet.clear();
				// System.out.println("Available: " + availableChoices);
				for (int q = 0; q < n; q++) {
					// Select the qth index from the available pool to use here:
					int randIndex = random.nextInt(n - q);
					// Find out what number this corresponds to, and write it in:
					Integer nextSelection = availableChoices.remove(randIndex);
					sets[s][q] = nextSelection.intValue();
				}
				// Track the chosen integers (in their selected order) to avoid duplicates
				for (int q = 0; q < n; q++) {
					// And track it in thisSet for hashing and 
					//  adding back to the pool
					thisSet.add(new Integer(sets[s][q]));
				}
				if (chosenSets.get(thisSet) == null) {
					// We haven't included this set yet:
					chosenSets.put(thisSet, new Integer(0));
					// System.out.println(" Chosen: " + thisSet);
					break;
				}
				// else we have already added this set,
				//  so we need to try again.
				// System.out.println(" Attempted: " + thisSet);
				// System.out.println(" Need to try again");
			}
		}
		return sets;
	}
	
	/**
	 * Generate all n! perturbations (assuming this doesn't blow our memory
	 * 
	 * @param n
	 * @return
	 */
	public int[][] generateAllDistinctPerturbations(int n) throws Exception {
		int numSets;
		try {
			numSets = MathsUtils.factorialCheckBounds(n);
		} catch (Exception e) {
			// n! blew Integer.MAX_INT
			throw new Exception("n! too large");
		}
		// allocate space for the distinct sets
		int[][] sets = new int[numSets][n];
		
		int[] workingSet = new int[n];
		Vector<Integer> availableChoices = new Vector<Integer>();
		for (int i = 0; i < n; i++) {
			availableChoices.add(new Integer(i));
		}
		addToDistinctPerturbations(sets, n, 0, workingSet, 0, availableChoices);
		return sets;
	}

	/**
	 * Using the workingSet which is filled up to (but not including)
	 *  fromIndex, add new distinct sets of nCp to the sets matrix,
	 *  from the setNumber index
	 * 
	 * @param sets
	 * @param n
	 * @param setNumber
	 * @param workingSet
	 * @param fromIndex
	 */
	protected int addToDistinctPerturbations(int[][] sets, int n, int setNumber, 
			int[] workingSet, int fromIndex, Vector<Integer> availableChoices) {
		if (fromIndex == n) {
			// The workingSet is ready to go, so copy it in
			// MatrixUtils.printArray(System.out, workingSet);
			System.arraycopy(workingSet, 0, sets[setNumber], 0, n);
			setNumber++;
		} else {
			// Iterate over a copy of the available choices, in case altering it during the loop
			//  causes problems.
			Vector<Integer> copyOfAvailableChoices = (Vector<Integer>) availableChoices.clone();
			// Add to the working set and pass it on:
			for (Integer nextInteger : copyOfAvailableChoices) {
				int nextInt = nextInteger.intValue();
				workingSet[fromIndex] = nextInt;
				// Remove this element as an available choice
				availableChoices.remove(nextInteger);
				// And keep filling out the array
				setNumber = addToDistinctPerturbations(sets, n, setNumber,
						workingSet, fromIndex + 1, availableChoices);
				// Put this integer back in as an available choice
				availableChoices.add(nextInteger);
			}
		}
		return setNumber;
	}

	/**
	 * Generate numberOfPerturbations perturbations of [0..n-1],
	 * which are not necessarily distinct.
	 * Could have double-ups even where the caller has asked for less
	 * than the number of distinct perturbations that exist.
	 * 
	 * @param n
	 * @param numberOfPerturbations
	 * @return an array of dimensions [numberOfPerturbations][n], with each row
	 *  being one perturbation of the elements
	 */
	public int[][] generateRandomPerturbations(int n, int numberOfPerturbations) {
		
		int[][] sets = new int[numberOfPerturbations][n];
		
		/* Manual implementation:
		for (int s = 0; s < numberOfPerturbations; s++) {
			// Generate a list of n random numbers:
			double[] randomList = generateRandomData(n);
			int[] sortedIndices = MatrixUtils.sortIndices(randomList);
			sets[s] = sortedIndices;
		}
		return sets;
		*/
		
		// Better implementation: using native classes, supplying
		// our Random object to ensure repeatability with the seed:
		
		// Use an array list because it gives RandomAccess to
		//  the Collections.shuffle method:
		ArrayList<Integer> list = new ArrayList<Integer>();
		for (int i = 0; i < n; i++) {
			list.add(i);
		}
		for (int s = 0; s < numberOfPerturbations; s++) {
			// Perform linear time shuffles (of what was already shuffled),
			//  Note: the shuffles are all equal likelihood
			Collections.shuffle(list, random);
			for (int j = 0; j < n; j++) {
				sets[s][j] = list.get(j);
			}
		}
		return sets;
	}

	public static void main(String[] args) throws Exception {
		// This code demonstrates that the Hashtable is hashing the 
		//  pointer rather than the array values -
		//  actually, it's hard to say, but if it is hashing the values then it's not doing .equals
		//  properly on the array because it isn't implemented.
		//  Should use a vector or an array object wrapper.
		java.util.Hashtable<int[],Integer> hashtable = new java.util.Hashtable<int[],Integer>();
		int[] array1 = {1,2,3,4,5};
		int[] array2 = {1,2,3,4,5};
		int[] array3 = {1,2,3,5,5};
		hashtable.put(array1, 0);
		hashtable.put(array3, 1);
		System.out.println(hashtable.get(array1));
		System.out.println(hashtable.get(array2));
		System.out.println(hashtable.get(array3));
		// This demonstrates is works on the values if we use vectors:
		Vector<Integer> vec1 = new Vector();
		Vector<Integer> vec2 = new Vector();
		Vector<Integer> vec3 = new Vector();
		for (int i = 0; i < array1.length; i++) {
			vec1.add(new Integer(array1[i]));
			vec2.add(new Integer(array2[i]));
			vec3.add(new Integer(array3[i]));
		}
		java.util.Hashtable<Vector<Integer>,Integer> hashtable2 = new java.util.Hashtable<Vector<Integer>,Integer>();
		hashtable2.put(vec1, 1);
		hashtable2.put(vec3, 3);
		System.out.println(hashtable2.get(vec1));
		System.out.println(hashtable2.get(vec2));
		System.out.println(hashtable2.get(vec3));
		
		RandomGenerator rg = new RandomGenerator();
		// And test out generating distinct random set:
		//MatrixUtils.printMatrix(System.out, rg.generateDistinctRandomSets(5, 2, 9));
		// Test out generating whole sets:
		MatrixUtils.printMatrix(System.out, rg.generateAllDistinctSets(5,3));

		System.out.println("Generating all distinct perturbations of 5:");
		MatrixUtils.printMatrix(System.out, rg.generateAllDistinctPerturbations(5));
		System.out.println("Generating 10 distinct perturbations of 4:");
		MatrixUtils.printMatrix(System.out, rg.generateDistinctRandomPerturbations(4, 10));
	}
	
	public class RandomPairs {
		public int n1, n2, p1, p2, N;
		public int[][] sets1;
		public int[][] sets2;
		
		public RandomPairs(int n1, int n2, int p1, int p2, int N) {
			this.n1 = n1;
			this.n2 = n2;
			this.p1 = p1;
			this.p2 = p2;
			this.N = N;
			sets1 = null;
			sets2 = null;
		}
	}
	
	/**
	 * Generate up to N <b>distinct<b/> pairs of p1 numbers from [0 .. n1-1]
	 *  and p2 numbers from [0 .. n2 - 1].
	 * 
	 * @param n1
	 * @param n2
	 * @param p1
	 * @param p2
	 * @param N
	 * @return 
	 */
	public RandomPairs generateDistinctPairsOfRandomSets(int n1, int n2, int p1, int p2,
			int N) {
		RandomPairs randPairs = new RandomPairs(n1, n2, p1, p2, N);
		int numOfPossibleSets1, numOfPossibleSets2;
		try {
			numOfPossibleSets1 = MathsUtils.numOfSets(n1, p1);
		} catch (Exception e) {
			// n1 choose p1 blew Integer.MAX_INT
			numOfPossibleSets1 = Integer.MAX_VALUE;
		}
		try {
			numOfPossibleSets2 = MathsUtils.numOfSets(n2, p2);
		} catch (Exception e) {
			// n2 choose p2 blew Integer.MAX_INT
			numOfPossibleSets2 = Integer.MAX_VALUE;
		}

		int[][] sets1 = generateDistinctRandomSets(n1, p1, N);
		int[][] sets2 = generateDistinctRandomSets(n2, p2, N);
		if ((numOfPossibleSets1 < N) || (numOfPossibleSets2 < N)) {
			// One pair does not have enough.
			// Use a long to avoid overflow
			long totalPossiblePairs = numOfPossibleSets1 * numOfPossibleSets2;
			if (totalPossiblePairs < N) {
				// We can return the product of the pairs
				randPairs.sets1 = new int[(int) totalPossiblePairs][];
				randPairs.sets2 = new int[(int) totalPossiblePairs][];
				int pairIndex = 0;
				for (int i1 = 0; i1 < sets1.length; i1++) {
					for (int i2 = 0; i2 < sets2.length; i2++) {
						randPairs.sets1[pairIndex] = sets1[i1];
						randPairs.sets2[pairIndex] = sets2[i2];
						pairIndex++;
					}
				}
			} else {
				// Need to randomly select pairs out of the ones we've already got here
				randPairs.sets1 = new int[N][];
				randPairs.sets2 = new int[N][];
				// Hashtable tracking the pairs we've already chosen
				HashSet<Vector<Integer>> alreadyChosen = new HashSet<Vector<Integer>>();
				for (int i = 0; i < N; i++) {
					// Select one of each of the distinct sets at random
					//  for the i-th pair
					for (;;) {
						// Select a candidate pair
						Vector<Integer> candidate = new Vector<Integer>();
						int candidate1 = random.nextInt(sets1.length);
						int candidate2 = random.nextInt(sets2.length);
						candidate.clear();
						candidate.add(candidate1);
						candidate.add(candidate2);
						if (!alreadyChosen.contains(candidate)) {
							// We're clear to add this candidate pair
							randPairs.sets1[i] = sets1[candidate1];
							randPairs.sets2[i] = sets2[candidate2];
							alreadyChosen.add(candidate);
							break;
						}
					}
				}
			}
		} else {
			// Both sets have enough pairs, so we can just use these without any repeats
			randPairs.sets1 = sets1;
			randPairs.sets2 = sets2;
		}
		return randPairs;
	}
}
