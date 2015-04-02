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

import java.util.Calendar;
import java.util.Collection;
import java.util.PriorityQueue;
import java.util.Vector;

import infodynamics.utils.KdTree.KdTreeNode;
import junit.framework.TestCase;

public class KdTreeTest extends TestCase {

	RandomGenerator rg = new RandomGenerator();
	
	public void testSmallConstruction() {
		// Testing with an example from 
		// http://en.wikipedia.org/wiki/K-d_tree
		
		double[][] data = { {2,3}, {5,4}, {9,6}, {4,7}, {8,1}, {7,2} };
		
		// KdTree kdTree = new KdTree( new int[] {2} , new double[][][] {data} );
		KdTree kdTree = new KdTree(data);
		//kdTree.print();
		validateAllPointsInTree(kdTree, data);
		validateStructure(kdTree, data);
	}

	public void testMultilayerConstruction() {
		
		double[][] data = { {1,15}, {2,14}, {3,13}, {4,12}, {5,11}, {6,10},
				{7,9}, {8,8}, {9,7}, {10,6}, {11,5}, {12,4}, {13,3},
				{14,2}, {15,1}};
		
		// KdTree kdTree = new KdTree( new int[] {2} , new double[][][] {data} );
		KdTree kdTree = new KdTree(data);
		//kdTree.print();
		validateAllPointsInTree(kdTree, data);
		validateStructure(kdTree, data);
	}
	
	
	public void testLargerConstruction() {
		int dimension = 4;
		double[][] data = rg.generateNormalData(100000, dimension, 0, 1);
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		// KdTree kdTree = new KdTree( new int[] {2} , new double[][][] {data} );
		KdTree kdTree = new KdTree(data);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Tree of %d points constructed in: %.3f sec\n",
				data.length, ((double) (endTimeTree - startTime)/1000.0));
		
		startTime = Calendar.getInstance().getTimeInMillis();
		validateAllPointsInTree(kdTree, data);
		validateStructure(kdTree, data);
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All points found and structure validated in: %.3f sec\n",
				((double) (endTimeValidate - startTime)/1000.0));
	}
	
	public void validateAllPointsInTree(KdTree kdtree, double[][] listOfPoints) {
		for (int t = 0; t < listOfPoints.length; t++) {
			KdTreeNode node = findNode(kdtree.rootNode, t, listOfPoints, 0);
			assertTrue(node != null);
		}
	}
	
	public KdTreeNode findNode(KdTreeNode startNode, int index, double[][] listOfPoints,
			int level) {
		if (startNode == null) {
			return null;
		}
		if (startNode.indexOfThisPoint == index) {
			return startNode;
		}
		if (listOfPoints[index][level % listOfPoints[index].length] <
				listOfPoints[startNode.indexOfThisPoint][level % listOfPoints[index].length]) {
			// search the left subtree
			return findNode(startNode.leftTree, index, listOfPoints, level + 1);
		} else {
			// search the right subtree
			return findNode(startNode.rightTree, index, listOfPoints, level + 1);
		}
	}
	
	public void validateStructure(KdTree kdTree, double[][] listOfPoints) {
		validateStructure(kdTree.rootNode, listOfPoints, 0, kdTree.totalDimensions,
				new Vector<Integer>(), new Vector<Double>(),
				new Vector<Integer>(), new Vector<Double>());
	}
	
	/**
	 * Validates that this node and its children meet the supplied
	 *  constraints, and also that all left children of this node
	 *  have corresponding value for the correct dimension strictly less
	 *  than this node, and all right children of this node have 
	 *  corresponding value for the correct dimension strictly greater than
	 *  this node.
	 * 
	 * @param startNode
	 * @param listOfPoints
	 * @param level
	 * @param totalDimensions
	 * @param indicesToBeLessThan
	 * @param lessThanValues
	 * @param indicesToBeGeq
	 * @param geqValues
	 */
	public void validateStructure(KdTreeNode startNode, double[][] listOfPoints,
			int level, int totalDimensions,
			Vector<Integer> indicesToBeLessThan, Vector<Double> lessThanValues,
			Vector<Integer> indicesToBeGeq, Vector<Double> geqValues) {
		
		if (startNode == null) {
			return;
		}
		
		int indexInValues = 0;
		for (Integer index : indicesToBeLessThan) {
			// Make sure the value for this node is less than
			//  the given value:
			assert(listOfPoints[startNode.indexOfThisPoint][index] <
					lessThanValues.get(indexInValues++));
		}
		indexInValues = 0;
		for (Integer index : indicesToBeGeq) {
			// Make sure the value for this node is >= 
			//  the given value:
			assert(listOfPoints[startNode.indexOfThisPoint][index] <
					geqValues.get(indexInValues++));
		}
		
		// Now check the children:
		// Left children require index level+1 % totalDimensions
		//  to be strictly less than that of this current point:
		if (startNode.leftTree != null) {
			indicesToBeLessThan.add((level+1) % totalDimensions);
			lessThanValues.add(
					listOfPoints[startNode.indexOfThisPoint][(level+1) % totalDimensions]);
			validateStructure(startNode.leftTree, listOfPoints, level+1, totalDimensions,
					indicesToBeLessThan, lessThanValues, indicesToBeGeq, geqValues);
			indicesToBeLessThan.remove(indicesToBeLessThan.lastElement());
			lessThanValues.remove(lessThanValues.lastElement());
		}
		// Right children require index level+1 % totalDimensions
		//  to be >= to that of this current point:
		if (startNode.rightTree != null) {
			indicesToBeGeq.add((level+1) % totalDimensions);
			geqValues.add(
					listOfPoints[startNode.indexOfThisPoint][(level+1) % totalDimensions]);
			validateStructure(startNode.rightTree, listOfPoints, level+1, totalDimensions,
					indicesToBeLessThan, lessThanValues, indicesToBeGeq, geqValues);
			indicesToBeGeq.remove(indicesToBeGeq.lastElement());
			geqValues.remove(geqValues.lastElement());
		}
	}

	public void testLargerConstructionWithDuplicates() {
		int dimension = 4;
		double[][] dataRaw = rg.generateNormalData(50000, dimension, 0, 1);
		double[][] data = new double[dataRaw.length * 2][];
		for (int t = 0; t < dataRaw.length; t++) {
			data[t] = dataRaw[t];
		}
		// Duplicate the first half into the second half here
		for (int t = dataRaw.length; t < 2*dataRaw.length; t++) {
			data[t] = dataRaw[t-dataRaw.length];
		}
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		// KdTree kdTree = new KdTree( new int[] {2} , new double[][][] {data} );
		KdTree kdTree = new KdTree(data);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Tree of %d points constructed in: %.3f sec\n",
				data.length, ((double) (endTimeTree - startTime)/1000.0));
		
		startTime = Calendar.getInstance().getTimeInMillis();
		validateAllPointsInTree(kdTree, data);
		validateStructure(kdTree, data);
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All points found and structure validated in: %.3f sec\n",
				((double) (endTimeValidate - startTime)/1000.0));
	}
	
	public void testFindNearestNeighbour() {
		int dimension = 4;
		double[][] data = rg.generateNormalData(1000, dimension, 0, 1);
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		KdTree kdTree = new KdTree(data);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Tree of %d points constructed in: %.3f sec\n",
				data.length, ((double) (endTimeTree - startTime)/1000.0));
		
		EuclideanUtils normCalculator = new EuclideanUtils(EuclideanUtils.NORM_MAX_NORM);
		startTime = Calendar.getInstance().getTimeInMillis();
		for (int t = 0; t < data.length; t++) {
			NeighbourNodeData nnData = kdTree.findNearestNeighbour(t);
			assertTrue(nnData != null);
			// Now find the nearest neighbour with a naive all-pairs comparison
			double currentMin = Double.POSITIVE_INFINITY;
			int currentMinIndex = 0;
			for (int t2 = 0; t2 < data.length; t2++) {
				double norm = Double.POSITIVE_INFINITY;
				if (t2 != t) {
					norm = normCalculator.norm(data[t], data[t2]);
				}
				if (norm < currentMin) {
					currentMin = norm;
					currentMinIndex = t2;
				}
			}
			if (currentMinIndex != nnData.sampleIndex) {
				System.out.printf("All pairs    : index %d, distance %.3f\n",
						currentMinIndex, currentMin);
				System.out.printf("kdTree search: index %d, distance %.3f\n",
						nnData.sampleIndex, nnData.distance);
			}
			assertEquals(currentMinIndex, nnData.sampleIndex);
		}		
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All nearest neighbours found in: %.3f sec\n",
				((double) (endTimeValidate - startTime)/1000.0));
	}

	public void testFindNearestNeighbourSeparateArrays() {
		int variables = 3;
		int dimensionsPerVariable = 3;
		int samples = 1000;
		double[][][] data = new double[variables][][];
		for (int v = 0; v < variables; v++) {
			data[v] = rg.generateNormalData(samples, dimensionsPerVariable, 0, 1);
		}
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		int[] dimensions = new int[variables];
		for (int v = 0; v < variables; v++) {
			dimensions[v] = dimensionsPerVariable;
		}
		KdTree kdTree = new KdTree(dimensions, data);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Tree of %d points constructed in: %.3f sec\n",
				samples, ((double) (endTimeTree - startTime)/1000.0));
		
		verifyNearestNeighbourForSeparateArrays(EuclideanUtils.NORM_MAX_NORM,
				data, samples, kdTree);
		// Must test with EuclideanUtils.NORM_EUCLIDEAN_SQUARED instead
		//  of EuclideanUtils.NORM_EUCLIDEAN if we want the min distances
		//  to match when tested inside verifyNearestNeighbourForSeparateArrays()
		verifyNearestNeighbourForSeparateArrays(EuclideanUtils.NORM_EUCLIDEAN_SQUARED,
				data, samples, kdTree);
	}
	
	private void verifyNearestNeighbourForSeparateArrays(int normType,
			double[][][] data, int samples, KdTree kdTree) {
		
		int variables = data.length;
		
		// Set up the given norm type:
		EuclideanUtils normCalculator = new EuclideanUtils(normType);
		kdTree.setNormType(normType);
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		for (int t = 0; t < samples; t++) {
			NeighbourNodeData nnData = kdTree.findNearestNeighbour(t);
			assertTrue(nnData != null);
			// Now find the nearest neighbour with a naive all-pairs comparison
			double currentMin = Double.POSITIVE_INFINITY;
			int currentMinIndex = 0;
			for (int t2 = 0; t2 < samples; t2++) {
				double norm = Double.POSITIVE_INFINITY;
				if (t2 != t) {
					double maxNorm = 0;
					for (int v = 0; v < variables; v++) {
						double normForThisVariable = normCalculator.norm(
								data[v][t], data[v][t2]);
						if (maxNorm < normForThisVariable) {
							maxNorm = normForThisVariable;
						}
					}
					norm = maxNorm;
				}
				if (norm < currentMin) {
					currentMin = norm;
					currentMinIndex = t2;
				}
			}
			if (currentMinIndex != nnData.sampleIndex) {
				System.out.printf("All pairs    : index %d, distance %.3f\n",
						currentMinIndex, currentMin);
				System.out.printf("kdTree search: index %d, distance %.3f\n",
						nnData.sampleIndex, nnData.distance);
			}
			assertEquals(currentMinIndex, nnData.sampleIndex);
			assertEquals(currentMin, nnData.distance);
		}		
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All nearest neighbours found in: %.3f sec\n",
				((double) (endTimeValidate - startTime)/1000.0));
	}

	public void testFindKNearestNeighbours() throws Exception {
		int dimension = 4;
		
		for (int K = 1; K < 5; K++) {
			double[][] data = rg.generateNormalData(1000, dimension, 0, 1);
			
			long startTime = Calendar.getInstance().getTimeInMillis();
			KdTree kdTree = new KdTree(data);
			long endTimeTree = Calendar.getInstance().getTimeInMillis();
			System.out.printf("Tree of %d points for %d NNs constructed in: %.3f sec\n",
					data.length, K, ((double) (endTimeTree - startTime)/1000.0));
			
			EuclideanUtils normCalculator = new EuclideanUtils(EuclideanUtils.NORM_MAX_NORM);
			startTime = Calendar.getInstance().getTimeInMillis();
			for (int t = 0; t < data.length; t++) {
				PriorityQueue<NeighbourNodeData> nnPQ =
						kdTree.findKNearestNeighbours(K, t);
				assertTrue(nnPQ.size() == K);
				// Also check the result is still correct if we search with
				//  zero size exclusion window:
				PriorityQueue<NeighbourNodeData> nnPQ_zeroExclWindow =
						kdTree.findKNearestNeighbours(K, t, 0);
				assertTrue(nnPQ_zeroExclWindow.size() == K);
				// Now find the K nearest neighbours with a naive all-pairs comparison
				double[][] distancesAndIndices = new double[data.length][2];
				for (int t2 = 0; t2 < data.length; t2++) {
					if (t2 != t) {
						distancesAndIndices[t2][0] = normCalculator.norm(data[t], data[t2]);
					} else {
						distancesAndIndices[t2][0] = Double.POSITIVE_INFINITY;
					}
					distancesAndIndices[t2][1] = t2;
				}
				int[] timeStepsOfKthMins =
						MatrixUtils.kMinIndices(distancesAndIndices, 0, K);
				for (int i = 0; i < K; i++) {
					// Check that the ith nearest neighbour matches for each method.
					// Note that these two method provide a different sorting order
					NeighbourNodeData nnData = nnPQ.poll();
					NeighbourNodeData nnData_zeroWindow = nnPQ_zeroExclWindow.poll();
					if (timeStepsOfKthMins[K - 1 - i] != nnData.sampleIndex) {
						// We have an error:
						System.out.printf("Erroneous match between indices %d (expected) " +
						 " and %d\n", timeStepsOfKthMins[K - 1 - i], nnData.sampleIndex);
					}
					assertEquals(timeStepsOfKthMins[K - 1 - i], nnData.sampleIndex);
					assertEquals(timeStepsOfKthMins[K - 1 - i], nnData_zeroWindow.sampleIndex);
				}
			}		
			long endTimeValidate = Calendar.getInstance().getTimeInMillis();
			System.out.printf("All %d nearest neighbours found in: %.3f sec\n",
					K, ((double) (endTimeValidate - startTime)/1000.0));
		}
	}

	public void testFindKNearestNeighboursWithExclusionWindow() throws Exception {
		int dimension = 4;
		int numSamples = 1000;
		int exclusionWindow = 100;
		
		for (int K = 1; K < 5; K++) {
			double[][] data = rg.generateNormalData(numSamples, dimension, 0, 1);
			
			long startTime = Calendar.getInstance().getTimeInMillis();
			KdTree kdTree = new KdTree(data);
			long endTimeTree = Calendar.getInstance().getTimeInMillis();
			System.out.printf("Tree of %d points for %d NNs constructed in: %.3f sec\n",
					data.length, K, ((double) (endTimeTree - startTime)/1000.0));
			
			EuclideanUtils normCalculator = new EuclideanUtils(EuclideanUtils.NORM_MAX_NORM);
			startTime = Calendar.getInstance().getTimeInMillis();
			for (int t = 0; t < data.length; t++) {
				PriorityQueue<NeighbourNodeData> nnPQ =
						kdTree.findKNearestNeighbours(K, t, exclusionWindow);
				assertTrue(nnPQ.size() == K);
				// Now find the K nearest neighbours with a naive all-pairs comparison
				double[][] distancesAndIndices = new double[data.length][2];
				for (int t2 = 0; t2 < data.length; t2++) {
					if (Math.abs(t2 - t) > exclusionWindow) {
						distancesAndIndices[t2][0] = normCalculator.norm(data[t], data[t2]);
					} else {
						distancesAndIndices[t2][0] = Double.POSITIVE_INFINITY;
					}
					distancesAndIndices[t2][1] = t2;
				}
				int[] timeStepsOfKthMins =
						MatrixUtils.kMinIndices(distancesAndIndices, 0, K);
				for (int i = 0; i < K; i++) {
					// Check that the ith nearest neighbour matches for each method.
					// Note that these two method provide a different sorting order
					NeighbourNodeData nnData = nnPQ.poll();
					if (timeStepsOfKthMins[K - 1 - i] != nnData.sampleIndex) {
						// We have an error:
						System.out.printf("Erroneous match between indices %d (expected) " +
						 " and %d\n", timeStepsOfKthMins[K - 1 - i], nnData.sampleIndex);
					}
					assertEquals(timeStepsOfKthMins[K - 1 - i], nnData.sampleIndex);
				}
			}		
			long endTimeValidate = Calendar.getInstance().getTimeInMillis();
			System.out.printf("All %d nearest neighbours found in: %.3f sec\n",
					K, ((double) (endTimeValidate - startTime)/1000.0));
		}
	}

	public void testFindKNearestNeighboursForSeparateArrays() throws Exception {
		int variables = 3;
		int dimensionsPerVariable = 3;
		int samples = 1000;
		int[] dimensions = new int[variables];
		for (int v = 0; v < variables; v++) {
			dimensions[v] = dimensionsPerVariable;
		}
		
		for (int K = 1; K < 5; K++) {
			double[][][] data = new double[variables][][];
			for (int v = 0; v < variables; v++) {
				data[v] = rg.generateNormalData(samples, dimensionsPerVariable, 0, 1);
			}
			
			long startTime = Calendar.getInstance().getTimeInMillis();
			KdTree kdTree = new KdTree(dimensions, data);
			long endTimeTree = Calendar.getInstance().getTimeInMillis();
			System.out.printf("Tree of %d points for %d NNs constructed in: %.3f sec\n",
					samples, K, ((double) (endTimeTree - startTime)/1000.0));
			
			verifyKNearestNeighboursForSeparateArrays(EuclideanUtils.NORM_MAX_NORM,
				data, K, samples, kdTree);
			// Must test with EuclideanUtils.NORM_EUCLIDEAN_SQUARED instead
			//  of EuclideanUtils.NORM_EUCLIDEAN if we want the min distances
			//  to match when tested inside verifyNearestNeighbourForSeparateArrays()
			verifyKNearestNeighboursForSeparateArrays(EuclideanUtils.NORM_EUCLIDEAN_SQUARED,
				data, K, samples, kdTree);
		}
	}

	private void verifyKNearestNeighboursForSeparateArrays(int normType,
			double[][][] data, int K, int samples, KdTree kdTree) throws Exception {
		
		int variables = data.length;
		
		// Set up the given norm type:
		EuclideanUtils normCalculator = new EuclideanUtils(normType);
		kdTree.setNormType(normType);

		long startTime = Calendar.getInstance().getTimeInMillis();
		for (int t = 0; t < samples; t++) {
			PriorityQueue<NeighbourNodeData> nnPQ =
					kdTree.findKNearestNeighbours(K, t);
			assertTrue(nnPQ.size() == K);
			// Now find the K nearest neighbours with a naive all-pairs comparison
			double[][] distancesAndIndices = new double[samples][2];
			for (int t2 = 0; t2 < samples; t2++) {
				if (t2 != t) {
					double maxNorm = 0;
					for (int v = 0; v < variables; v++) {
						double normForThisVariable = normCalculator.norm(
								data[v][t], data[v][t2]);
						if (maxNorm < normForThisVariable) {
							maxNorm = normForThisVariable;
						}
					}
					distancesAndIndices[t2][0] = maxNorm;
				} else {
					distancesAndIndices[t2][0] = Double.POSITIVE_INFINITY;
				}
				distancesAndIndices[t2][1] = t2;
			}
			int[] timeStepsOfKthMins =
					MatrixUtils.kMinIndices(distancesAndIndices, 0, K);
			for (int i = 0; i < K; i++) {
				// Check that the ith nearest neighbour matches for each method.
				// Note that these two method provide a different sorting order
				NeighbourNodeData nnData = nnPQ.poll();
				if (timeStepsOfKthMins[K - 1 - i] != nnData.sampleIndex) {
					// We have an error:
					System.out.printf("Erroneous match between indices %d (expected) " +
					 " and %d\n", timeStepsOfKthMins[K - 1 - i], nnData.sampleIndex);
				}
				assertEquals(timeStepsOfKthMins[K - 1 - i], nnData.sampleIndex);
			}
		}		
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All %d nearest neighbours found in: %.3f sec\n",
				K, ((double) (endTimeValidate - startTime)/1000.0));
	}

	public void testCountNeighboursWithinRSeparateArrays() {
		int variables = 3;
		int dimensionsPerVariable = 3;
		int samples = 1000;
		double[][][] data = new double[variables][][];
		for (int v = 0; v < variables; v++) {
			data[v] = rg.generateNormalData(samples, dimensionsPerVariable, 0, 1);
		}
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		int[] dimensions = new int[variables];
		for (int v = 0; v < variables; v++) {
			dimensions[v] = dimensionsPerVariable;
		}
		KdTree kdTree = new KdTree(dimensions, data);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Tree of %d points constructed in: %.3f sec\n",
				samples, ((double) (endTimeTree - startTime)/1000.0));
		
		double[] rs = {0.2, 0.6, 0.8, 1.0};
		for (int i = 0; i < rs.length; i++) {
			verifyCountNeighboursWithinRForSeparateArrays(EuclideanUtils.NORM_MAX_NORM,
					data, samples, kdTree, rs[i], true);
			verifyCountNeighboursWithinRForSeparateArrays(EuclideanUtils.NORM_MAX_NORM,
					data, samples, kdTree, rs[i], false);
			// And test the array return variant also:
			verifyCountNeighboursWithinRForSeparateArraysArrayArgs(EuclideanUtils.NORM_MAX_NORM,
					data, samples, kdTree, rs[i], true);
			verifyCountNeighboursWithinRForSeparateArraysArrayArgs(EuclideanUtils.NORM_MAX_NORM,
					data, samples, kdTree, rs[i], false);
			// Must test with EuclideanUtils.NORM_EUCLIDEAN_SQUARED instead
			//  of EuclideanUtils.NORM_EUCLIDEAN if we want the min distances
			//  to match when tested inside verifyNearestNeighbourForSeparateArrays()
			verifyCountNeighboursWithinRForSeparateArrays(EuclideanUtils.NORM_EUCLIDEAN_SQUARED,
					data, samples, kdTree, rs[i], true);
			verifyCountNeighboursWithinRForSeparateArrays(EuclideanUtils.NORM_EUCLIDEAN_SQUARED,
					data, samples, kdTree, rs[i], false);
			verifyCountNeighboursWithinRForSeparateArraysArrayArgs(EuclideanUtils.NORM_EUCLIDEAN_SQUARED,
					data, samples, kdTree, rs[i], true);
			verifyCountNeighboursWithinRForSeparateArraysArrayArgs(EuclideanUtils.NORM_EUCLIDEAN_SQUARED,
					data, samples, kdTree, rs[i], false);
		}
	}

	private void verifyCountNeighboursWithinRForSeparateArrays(int normType,
			double[][][] data, int samples, KdTree kdTree, double r,
			boolean allowEqualToR) {
		
		int variables = data.length;
		
		// Set up the given norm type:
		EuclideanUtils normCalculator = new EuclideanUtils(normType);
		kdTree.setNormType(normType);
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		int totalCount = 0;
		int[] counts = new int[samples];
		for (int t = 0; t < samples; t++) {
			// int count = kdTree.countPointsWithinR(t, r, allowEqualToR);
			Collection<NeighbourNodeData> pointsWithinR =
					kdTree.findPointsWithinR(t, r, allowEqualToR);
			int count = pointsWithinR.size();
			assertTrue(count >= 0);
			counts[t] = count;
			totalCount += count;
		}
		long endTimeCount = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All neighbours within %.3f (average of %.3f) found in: %.3f sec\n",
				r, (double) totalCount / (double) samples,
				((double) (endTimeCount - startTime)/1000.0));
		// Now find the neighbour count with a naive all-pairs comparison
		for (int t = 0; t < samples; t++) {
			int naiveCount = 0;
			for (int t2 = 0; t2 < samples; t2++) {
				double norm = Double.POSITIVE_INFINITY;
				if (t2 != t) {
					double maxNorm = 0;
					for (int v = 0; v < variables; v++) {
						double normForThisVariable = normCalculator.norm(
								data[v][t], data[v][t2]);
						if (maxNorm < normForThisVariable) {
							maxNorm = normForThisVariable;
						}
					}
					norm = maxNorm;
				}
				if ((!allowEqualToR && (norm < r)) ||
					( allowEqualToR && (norm <= r))) {
					naiveCount++;
				}
			}
			if (naiveCount != counts[t]) {
				System.out.printf("All pairs    : count %d\n", naiveCount);
				System.out.printf("kdTree search: count %d\n", counts[t]);
			}
			assertEquals(naiveCount, counts[t]);
		}
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All neighbours within %.3f (average of %.3f) validated in: %.3f sec\n",
				r, (double) totalCount / (double) samples,
				((double) (endTimeValidate - endTimeCount)/1000.0));
	}

	private void verifyCountNeighboursWithinRForSeparateArraysArrayArgs(int normType,
			double[][][] data, int samples, KdTree kdTree, double r,
			boolean allowEqualToR) {
		
		int variables = data.length;
		
		// Set up the given norm type:
		EuclideanUtils normCalculator = new EuclideanUtils(normType);
		kdTree.setNormType(normType);
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		int totalCount = 0;
		int[] counts = new int[samples];
		boolean[] withinR = new boolean[samples];
		int[] indicesWithinR = new int[samples];
		for (int t = 0; t < samples; t++) {
			// int count = kdTree.countPointsWithinR(t, r, allowEqualToR);
			kdTree.findPointsWithinR(t, r, allowEqualToR, withinR, indicesWithinR);
			int count = 0;
			// Run through the list of points returned as within r and check them
			for (int t2 = 0; indicesWithinR[t2] != -1; t2++) {
				count++;
				assertTrue(withinR[indicesWithinR[t2]]);
				// Reset this marker
				withinR[indicesWithinR[t2]] = false;
				// Check that it's really within r
				double maxNorm = 0;
				for (int v = 0; v < variables; v++) {
					double normForThisVariable = normCalculator.norm(
							data[v][t], data[v][indicesWithinR[t2]]);
					if (maxNorm < normForThisVariable) {
						maxNorm = normForThisVariable;
					}
				}
				
				assertTrue((!allowEqualToR && (maxNorm < r)) ||
						( allowEqualToR && (maxNorm <= r)));
			}
			// Check that there were no other points marked as within r:
			int count2 = 0;
			for (int t2 = 0; t2 < samples; t2++) {
				if (withinR[t2]) {
					count2++;
				}
			}
			assertEquals(0, count2); // We set all the ones that should have been true to false already
			// Now reset the boolean array
			counts[t] = count;
			totalCount += count;
		}
		long endTimeCount = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All neighbours within %.3f (average of %.3f) found in: %.3f sec\n",
				r, (double) totalCount / (double) samples,
				((double) (endTimeCount - startTime)/1000.0));
		// Now find the neighbour count with a naive all-pairs comparison
		for (int t = 0; t < samples; t++) {
			int naiveCount = 0;
			for (int t2 = 0; t2 < samples; t2++) {
				double norm = Double.POSITIVE_INFINITY;
				if (t2 != t) {
					double maxNorm = 0;
					for (int v = 0; v < variables; v++) {
						double normForThisVariable = normCalculator.norm(
								data[v][t], data[v][t2]);
						if (maxNorm < normForThisVariable) {
							maxNorm = normForThisVariable;
						}
					}
					norm = maxNorm;
				}
				if ((!allowEqualToR && (norm < r)) ||
					( allowEqualToR && (norm <= r))) {
					naiveCount++;
				}
			}
			if (naiveCount != counts[t]) {
				System.out.printf("All pairs    : count %d\n", naiveCount);
				System.out.printf("kdTree search: count %d\n", counts[t]);
			}
			assertEquals(naiveCount, counts[t]);
		}
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All neighbours within %.3f (average of %.3f) validated in: %.3f sec\n",
				r, (double) totalCount / (double) samples,
				((double) (endTimeValidate - endTimeCount)/1000.0));
	}

	public void testCountNeighboursWithinRExclusionWindow() {
		int dimensionsPerVariable = 3;
		int samples = 1000;
		double[][] data = rg.generateNormalData(samples, dimensionsPerVariable, 0, 1);
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		KdTree kdTree = new KdTree(data);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Tree of %d points constructed in: %.3f sec\n",
				samples, ((double) (endTimeTree - startTime)/1000.0));
		
		double[] rs = {0.2, 0.6, 0.8, 1.0};
		for (int i = 0; i < rs.length; i++) {
			verifyCountNeighboursWithinRExclusionWindow(EuclideanUtils.NORM_MAX_NORM,
					data, samples, kdTree, rs[i], true);
			verifyCountNeighboursWithinRExclusionWindow(EuclideanUtils.NORM_MAX_NORM,
					data, samples, kdTree, rs[i], false);
			// and check using the array methods as well
			verifyCountNeighboursWithinRExclusionWindowArrayArgs(EuclideanUtils.NORM_MAX_NORM,
					data, samples, kdTree, rs[i], true);
			verifyCountNeighboursWithinRExclusionWindowArrayArgs(EuclideanUtils.NORM_MAX_NORM,
					data, samples, kdTree, rs[i], false);
		}
	}

	private void verifyCountNeighboursWithinRExclusionWindow(int normType,
			double[][] data, int samples, KdTree kdTree, double r,
			boolean allowEqualToR) {
		
		int exclWindow = 100;
		
		// Set up the given norm type:
		EuclideanUtils normCalculator = new EuclideanUtils(normType);
		kdTree.setNormType(normType);
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		int totalCount = 0;
		int[] counts = new int[samples];
		for (int t = 0; t < samples; t++) {
			int count =
					kdTree.countPointsWithinR(t, r, exclWindow, allowEqualToR);
			assertTrue(count >= 0);
			counts[t] = count;
			totalCount += count;
		}
		long endTimeCount = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All neighbours within %.3f (average of %.3f) found in: %.3f sec\n",
				r, (double) totalCount / (double) samples,
				((double) (endTimeCount - startTime)/1000.0));
		// Now find the neighbour count with a naive all-pairs comparison
		for (int t = 0; t < samples; t++) {
			int naiveCount = 0;
			for (int t2 = 0; t2 < samples; t2++) {
				double norm = Double.POSITIVE_INFINITY;
				if (Math.abs(t2 - t) > exclWindow) {
					norm = normCalculator.norm(
							data[t], data[t2]);
				}
				if ((!allowEqualToR && (norm < r)) ||
					( allowEqualToR && (norm <= r))) {
					naiveCount++;
				}
			}
			if (naiveCount != counts[t]) {
				// Add some debug prints:
				System.out.printf("All pairs    : count %d\n", naiveCount);
				System.out.printf("kdTree search: count %d\n", counts[t]);
				Collection<NeighbourNodeData> pointsWithinR =
						kdTree.findPointsWithinR(t, r, allowEqualToR);
				for (NeighbourNodeData nnData : pointsWithinR) {
					System.out.printf("kd: t=%d: time diff %d, norm %.3f\n",
							t, Math.abs(nnData.sampleIndex - t),
							normCalculator.norm(
									data[t], data[nnData.sampleIndex]));
				}
				// and compare to all pairs:
				for (int t2 = 0; t2 < samples; t2++) {
					double norm = Double.POSITIVE_INFINITY;
					if (Math.abs(t2 - t) > exclWindow) {
						norm = normCalculator.norm(
								data[t], data[t2]);
					}
					if ((!allowEqualToR && (norm < r)) ||
						( allowEqualToR && (norm <= r))) {
						System.out.printf("naive: t=%d: time diff %d, norm %.3f\n",
								t, Math.abs(t2 - t), norm);
					}
				}
			}
			assertEquals(naiveCount, counts[t]);
		}
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All neighbours within %.3f (average of %.3f) validated in: %.3f sec\n",
				r, (double) totalCount / (double) samples,
				((double) (endTimeValidate - endTimeCount)/1000.0));
	}

	private void verifyCountNeighboursWithinRExclusionWindowArrayArgs(int normType,
			double[][] data, int samples, KdTree kdTree, double r,
			boolean allowEqualToR) {
		
		int exclWindow = 100;
		
		// Set up the given norm type:
		EuclideanUtils normCalculator = new EuclideanUtils(normType);
		kdTree.setNormType(normType);
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		int totalCount = 0;
		int[] counts = new int[samples];
		boolean[] withinR = new boolean[samples];
		int[] indicesWithinR = new int[samples];
		for (int t = 0; t < samples; t++) {
			// int count = kdTree.countPointsWithinR(t, r, allowEqualToR);
			kdTree.findPointsWithinR(t, r, exclWindow, allowEqualToR,
					withinR, indicesWithinR);
			int count = 0;
			// Run through the list of points returned as within r and check them
			for (int t2 = 0; indicesWithinR[t2] != -1; t2++) {
				count++;
				assertTrue(withinR[indicesWithinR[t2]]);
				// Reset this marker
				withinR[indicesWithinR[t2]] = false;
				// Check that it's really within r
				double maxNorm = normCalculator.norm(
							data[t], data[indicesWithinR[t2]]);
				assertTrue((!allowEqualToR && (maxNorm < r)) ||
						( allowEqualToR && (maxNorm <= r)));
			}
			// Check that there were no other points marked as within r:
			int count2 = 0;
			for (int t2 = 0; t2 < samples; t2++) {
				if (withinR[t2]) {
					count2++;
				}
			}
			assertEquals(0, count2); // We set all the ones that should have been true to false already
			// Now reset the boolean array
			counts[t] = count;
			totalCount += count;
		}
		long endTimeCount = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All neighbours within %.3f (average of %.3f) found in: %.3f sec\n",
				r, (double) totalCount / (double) samples,
				((double) (endTimeCount - startTime)/1000.0));
		// Now find the neighbour count with a naive all-pairs comparison
		for (int t = 0; t < samples; t++) {
			int naiveCount = 0;
			for (int t2 = 0; t2 < samples; t2++) {
				double norm = Double.POSITIVE_INFINITY;
				if (Math.abs(t2 - t) > exclWindow) {
					norm = normCalculator.norm(
							data[t], data[t2]);
				}
				if ((!allowEqualToR && (norm < r)) ||
					( allowEqualToR && (norm <= r))) {
					naiveCount++;
				}
			}
			if (naiveCount != counts[t]) {
				// Add some debug prints:
				System.out.printf("All pairs    : count %d\n", naiveCount);
				System.out.printf("kdTree search: count %d\n", counts[t]);
				Collection<NeighbourNodeData> pointsWithinR =
						kdTree.findPointsWithinR(t, r, allowEqualToR);
				for (NeighbourNodeData nnData : pointsWithinR) {
					System.out.printf("kd: t=%d: time diff %d, norm %.3f\n",
							t, Math.abs(nnData.sampleIndex - t),
							normCalculator.norm(
									data[t], data[nnData.sampleIndex]));
				}
				// and compare to all pairs:
				for (int t2 = 0; t2 < samples; t2++) {
					double norm = Double.POSITIVE_INFINITY;
					if (Math.abs(t2 - t) > exclWindow) {
						norm = normCalculator.norm(
								data[t], data[t2]);
					}
					if ((!allowEqualToR && (norm < r)) ||
						( allowEqualToR && (norm <= r))) {
						System.out.printf("naive: t=%d: time diff %d, norm %.3f\n",
								t, Math.abs(t2 - t), norm);
					}
				}
			}
			assertEquals(naiveCount, counts[t]);
		}
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All neighbours within %.3f (average of %.3f) validated in: %.3f sec\n",
				r, (double) totalCount / (double) samples,
				((double) (endTimeValidate - endTimeCount)/1000.0));
	}

	public void testCountNeighboursWithinRSeparateArraysWithDuplicates() throws Exception {
		int variables = 3;
		int dimensionsPerVariable = 3;
		int samples = 2000;
		double[][][] data = new double[variables][samples][dimensionsPerVariable];
		for (int v = 0; v < variables; v++) {
			double[][] rawData = rg.generateNormalData(samples/2, dimensionsPerVariable, 0, 1);
			// Now add the duplicates in
			for (int t = 0; t < samples/2; t++) {
				for (int c = 0; c < dimensionsPerVariable; c++) {
					data[v][t][c] = rawData[t][c];
					data[v][t+samples/2][c] = rawData[t][c];
				}
			}
		}
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		int[] dimensions = new int[variables];
		for (int v = 0; v < variables; v++) {
			dimensions[v] = dimensionsPerVariable;
		}
		KdTree kdTree = new KdTree(dimensions, data);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Tree of %d points constructed in: %.3f sec\n",
				samples, ((double) (endTimeTree - startTime)/1000.0));

		for (int K = 2; K < 5; K++) {
			startTime = Calendar.getInstance().getTimeInMillis();			
			for (int t = 0; t < data.length; t++) {
				PriorityQueue<NeighbourNodeData> nnPQ =
						kdTree.findKNearestNeighbours(K, t);
				assertTrue(nnPQ.size() == K);
				// What is the maximum distance?
				double maxDistance = nnPQ.peek().distance;
				// Now check how many points we can count within or at maxDistance:
				int countWithinOrOn = kdTree.countPointsWithinOrOnR(t, maxDistance);
				int countStrictlyWithin = kdTree.countPointsStrictlyWithinR(t, maxDistance);
				
				// The following assumes no duplicates apart from those we designed
				//  (this should be true), and that no sets of points
				//  (apart from the desiged duplicates) are at the same norm
				//  (again, it is highly unlikely that this is not true).
				//
				// With our design of duplicates of all points, we should always
				//  have an odd number of points strictly within the box
				//  (given the duplicate of the point itself, then pairs of
				//  points within), and an even number within or on (adding in
				//  the pair of points on the boundary).
				if (K % 2 == 0) {
					// K is even -- expect K - 1 within and K + 1 on
					assertEquals(K-1, countStrictlyWithin);
					assertEquals(K+1, countWithinOrOn);
				} else {
					// K is odd -- expect K - 2 within and K on
					assertEquals(K-2, countStrictlyWithin);
					assertEquals(K, countWithinOrOn);
				}
			}
		}
	}

	public void testCountNeighboursWithinRsSeparateArrays() {
		int variables = 3;
		int dimensionsPerVariable = 3;
		int samples = 1000;
		double[][][] data = new double[variables][][];
		for (int v = 0; v < variables; v++) {
			data[v] = rg.generateNormalData(samples, dimensionsPerVariable, 0, 1);
		}
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		int[] dimensions = new int[variables];
		for (int v = 0; v < variables; v++) {
			dimensions[v] = dimensionsPerVariable;
		}
		KdTree kdTree = new KdTree(dimensions, data);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Tree of %d points constructed in: %.3f sec\n",
				samples, ((double) (endTimeTree - startTime)/1000.0));
		
		double[] r1s = {0.2, 0.6, 0.8, 1.0};
		double[] r2s = {0.25, 0.5, 0.6, 0.5};
		double[] r3s = {0.5, 0.15, 0.2, 0.8};
		for (int i = 0; i < r1s.length; i++) {
			verifyCountNeighboursWithinRForSeparateArrays(EuclideanUtils.NORM_MAX_NORM,
					data, samples, kdTree, new double[] {r1s[i], r2s[i], r3s[i]}, true);
			verifyCountNeighboursWithinRForSeparateArrays(EuclideanUtils.NORM_MAX_NORM,
					data, samples, kdTree, new double[] {r1s[i], r2s[i], r3s[i]}, false);
			// Must test with EuclideanUtils.NORM_EUCLIDEAN_SQUARED instead
			//  of EuclideanUtils.NORM_EUCLIDEAN if we want the min distances
			//  to match when tested inside verifyNearestNeighbourForSeparateArrays()
			verifyCountNeighboursWithinRForSeparateArrays(EuclideanUtils.NORM_EUCLIDEAN_SQUARED,
					data, samples, kdTree, new double[] {r1s[i], r2s[i], r3s[i]}, true);
			verifyCountNeighboursWithinRForSeparateArrays(EuclideanUtils.NORM_EUCLIDEAN_SQUARED,
					data, samples, kdTree, new double[] {r1s[i], r2s[i], r3s[i]}, false);
		}
	}

	private void verifyCountNeighboursWithinRForSeparateArrays(int normType,
			double[][][] data, int samples, KdTree kdTree, double[] rs,
			boolean allowEqualToR) {
		
		int variables = data.length;
		
		// Set up the given norm type:
		EuclideanUtils normCalculator = new EuclideanUtils(normType);
		kdTree.setNormType(normType);
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		int totalCount = 0;
		for (int t = 0; t < samples; t++) {
			int count = kdTree.countPointsWithinRs(t, rs, allowEqualToR);
			assertTrue(count >= 0);
			totalCount += count;
			// Now find the neighbour count with a naive all-pairs comparison
			int naiveCount = 0;
			for (int t2 = 0; t2 < samples; t2++) {
				if (t2 == t) {
					continue;
				}
				boolean withinBounds = true;
				for (int v = 0; v < variables; v++) {
					double normForThisVariable = normCalculator.norm(
							data[v][t], data[v][t2]);
					if ((allowEqualToR && (normForThisVariable > rs[v])) ||
						 (!allowEqualToR && (normForThisVariable >= rs[v]))) {
						withinBounds = false;
						break;
					}
				}
				if (withinBounds) {
					naiveCount++;
				}
			}
			if (naiveCount != count) {
				System.out.printf("All pairs    : count %d\n", naiveCount);
				System.out.printf("kdTree search: count %d\n", count);
			}
			assertEquals(naiveCount, count);
		}		
		long endTimeValidate = Calendar.getInstance().getTimeInMillis();
		System.out.printf("All neighbours within bounds (average of %.3f) found in: %.3f sec\n",
				(double) totalCount / (double) samples,
				((double) (endTimeValidate - startTime)/1000.0));
	}

	public void testAdditionalCriteria() throws Exception {
		validateAdditionalCriteria(3, 3, true);
		validateAdditionalCriteria(3, 3, false);
		// The following will include a little unit testing
		//  of UnivariateNearestNeighbourSearcher
		//  for the data0 in both cases, and the remaining data in 
		//  the last:
		validateAdditionalCriteria(3, 1, false);
		validateAdditionalCriteria(2, 1, false);
	}
	
	public void validateAdditionalCriteria(int variables,
			int dimensionsPerVariable, boolean allowEqualsR) throws Exception {
		int samples = 1000;
		double[][][] data = new double[variables][][];
		double[][][] dataAfter1 = new double[variables-1][][];
		for (int v = 0; v < variables; v++) {
			data[v] = rg.generateNormalData(samples, dimensionsPerVariable, 0, 1);
			if (v > 0) {
				dataAfter1[v-1] = data[v];
			}
		}
		double[][] data1 = data[0];

		long startTime = Calendar.getInstance().getTimeInMillis();
		int[] dimensions = new int[variables];
		for (int v = 0; v < variables; v++) {
			dimensions[v] = dimensionsPerVariable;
		}
		int[] dimensionsAfter1 = new int[variables-1];
		for (int v = 1; v < variables; v++) {
			dimensionsAfter1[v-1] = dimensionsPerVariable;
		}
		KdTree kdTree = new KdTree(dimensions, data);
		NearestNeighbourSearcher nnSearcherExceptVar1 =
				NearestNeighbourSearcher.create(dimensionsAfter1, dataAfter1);
		NearestNeighbourSearcher nnSearcherVar1 =
				NearestNeighbourSearcher.create(data1);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Additional criteria test: Tree of %d points constructed in: %.3f sec\n",
				samples, ((double) (endTimeTree - startTime)/1000.0));

		double[] r1 = {1.0};
		boolean[] isWithinR = new boolean[samples];
		int[] indicesWithinR = new int[samples+1];
		for (int i = 0; i < r1.length; i++) {
			for (int t = 0; t < samples; t++) {
				// Establish our ground truth:
				int count = kdTree.countPointsWithinR(t, r1[i], allowEqualsR);
				// Alternative if we need to debug:
				// Collection<NeighbourNodeData> coll = kdTree.findPointsWithinR(t, r1[i], allowEqualsR);
				// int count = coll.size();
				
				// Now search in variable 1 first, then use that 
				//  to help the others:
				nnSearcherVar1.findPointsWithinR(t, r1[i], allowEqualsR,
						isWithinR, indicesWithinR);
				int countUsingVar1 = kdTree.countPointsWithinR(t, r1[i],
						allowEqualsR, 0, isWithinR);
				assertEquals(count, countUsingVar1);
				
				// And then search the kdTree which does not have variable 1:
				int countAdditionalCriteria =
						nnSearcherExceptVar1.countPointsWithinR(t, r1[i],
								allowEqualsR, isWithinR);
				assertEquals(count, countAdditionalCriteria);
				
				// Finally, reset our boolean array:
				for (int t2 = 0; indicesWithinR[t2] != -1; t2++) {
					isWithinR[indicesWithinR[t2]] = false;
				}
			}
		}
	}
	
	public void testPreviouslyTestedVariableMethodMultipleRs() throws Exception {
		
		for (int vars = 2; vars <= 3; vars++) {
			for (int dimsPerVar = 1; dimsPerVar <= 3; dimsPerVar++) {
				for (int varToPretest = 0; varToPretest < vars; varToPretest++) {
					validatePreviouslyTestedVariableMethodMultipleRs(
							vars, dimsPerVar, varToPretest, true);
					validatePreviouslyTestedVariableMethodMultipleRs(
							vars, dimsPerVar, varToPretest, false);
				}
			}
		}
	}
	
	public void validatePreviouslyTestedVariableMethodMultipleRs(int variables,
			int dimensionsPerVariable, int variableToPretest, boolean allowEqualsR) throws Exception {
		int samples = 1000;
		double[][][] data = new double[variables][][];
		for (int v = 0; v < variables; v++) {
			data[v] = rg.generateNormalData(samples, dimensionsPerVariable, 0, 1);
		}

		long startTime = Calendar.getInstance().getTimeInMillis();
		int[] dimensions = new int[variables];
		for (int v = 0; v < variables; v++) {
			dimensions[v] = dimensionsPerVariable;
		}
		KdTree kdTree = new KdTree(dimensions, data);
		NearestNeighbourSearcher nnSearcherVarToPretest =
				NearestNeighbourSearcher.create(data[variableToPretest]);
		long endTimeTree = Calendar.getInstance().getTimeInMillis();
		System.out.printf("Previously Tested Variable Method Multiple Rs test: Tree of %d points constructed in: %.3f sec\n",
				samples, ((double) (endTimeTree - startTime)/1000.0));

		double[] mainRadii = {1.0};
		boolean[] isWithinR = new boolean[samples];
		int[] indicesWithinR = new int[samples+1];
		for (int i = 0; i < mainRadii.length; i++) {
			// Construct a set of radii here:
			double[] r = new double[variables];
			for (int v = 0; v < variables; v++) {
				r[v] = mainRadii[i] * Math.pow(0.7, v);
			}
			// Now test each sample point:
			for (int t = 0; t < samples; t++) {
				// Establish our ground truth:
				int count = kdTree.countPointsWithinRs(t, r, allowEqualsR);
				
				// Now search in the pretest first, then use that 
				//  to help the others:
				nnSearcherVarToPretest.findPointsWithinR(t, r[variableToPretest], 0,
						allowEqualsR, isWithinR, indicesWithinR);
				int countUsingPretestVar = kdTree.countPointsWithinRs(t, r,
						allowEqualsR, variableToPretest, isWithinR);
				assertEquals(count, countUsingPretestVar);
				
				// Finally, reset our boolean array:
				for (int t2 = 0; indicesWithinR[t2] != -1; t2++) {
					isWithinR[indicesWithinR[t2]] = false;
				}
			}
		}
	}

}
