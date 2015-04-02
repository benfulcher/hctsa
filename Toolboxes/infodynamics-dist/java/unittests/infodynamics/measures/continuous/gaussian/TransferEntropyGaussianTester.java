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

package infodynamics.measures.continuous.gaussian;

import java.util.Vector;

import infodynamics.measures.continuous.TransferEntropyAbstractTester;
import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

public class TransferEntropyGaussianTester extends TransferEntropyAbstractTester {

	public void testLocalsAverageCorrectly() throws Exception {
		TransferEntropyCalculatorGaussian teCalc = new TransferEntropyCalculatorGaussian();
		super.testLocalsAverageCorrectly(teCalc, 100, 2);
	}

	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		TransferEntropyCalculatorGaussian teCalc = new TransferEntropyCalculatorGaussian();
		super.testComputeSignificanceDoesntAlterAverage(teCalc, 100, 2);
	}

	/**
	 * Compare to Granger causality values for 2coupledRandomCols-1.txt
	 *  as computed by modifications to computeGranger.m from 
	 *  ChaLearn Connectomics Challenge Sample Code
 	 *   (http://connectomics.chalearn.org).
 	 * The modifications returned log(RSS0/RSS1) in that code to match
 	 *  Barnett and Seth's Granger definition which is 2 x TE. 
	 * 
	 * @throws Exception
	 */
	public void testGConTestData() throws Exception {
		System.out.println("Testing linear-Gaussian TE against values from " +
				"(modified) ChaLearn Connectomics Challenge Sample Code");
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/2coupledRandomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		TransferEntropyCalculatorGaussian teCalc = new TransferEntropyCalculatorGaussian();
		
		// Check TE(column 0 -> column 1) against Granger causality result (k,l=1):
		teCalc.initialise(1);
		teCalc.setObservations(MatrixUtils.selectColumn(data, 0),
				MatrixUtils.selectColumn(data, 1));
		double te = teCalc.computeAverageLocalOfObservations();
		assertEquals(0.71693/2, te, 0.0001);

		// Check TE(column 1 -> column 0) against Granger causality result  (k,l=1):
		teCalc.initialise(1);
		teCalc.setObservations(MatrixUtils.selectColumn(data, 1),
				MatrixUtils.selectColumn(data, 0));
		te = teCalc.computeAverageLocalOfObservations();
		assertEquals(0.01702/2, te, 0.0001);

		// Check TE(column 0 -> column 1) against Granger causality result (k,l=2):
		teCalc.initialise(2, 1, 2, 1, 1);
		teCalc.setObservations(MatrixUtils.selectColumn(data, 0),
				MatrixUtils.selectColumn(data, 1));
		te = teCalc.computeAverageLocalOfObservations();
		assertEquals(0.77806/2, te, 0.0001);

		// Check TE(column 1 -> column 0) against Granger causality result (k,l=2):
		teCalc.initialise(2, 1, 2, 1, 1);
		teCalc.setObservations(MatrixUtils.selectColumn(data, 1),
				MatrixUtils.selectColumn(data, 0));
		te = teCalc.computeAverageLocalOfObservations();
		assertEquals(0.02407/2, te, 0.0001);

		// Check TE(column 0 -> column 1) against Granger causality result (k,l=3):
		teCalc.initialise(3, 1, 3, 1, 1);
		teCalc.setObservations(MatrixUtils.selectColumn(data, 0),
				MatrixUtils.selectColumn(data, 1));
		te = teCalc.computeAverageLocalOfObservations();
		assertEquals(0.79468/2, te, 0.0001);

		// Check TE(column 1 -> column 0) against Granger causality result (k,l=3):
		teCalc.initialise(3, 1, 3, 1, 1);
		teCalc.setObservations(MatrixUtils.selectColumn(data, 1),
				MatrixUtils.selectColumn(data, 0));
		te = teCalc.computeAverageLocalOfObservations();
		assertEquals(0.02180/2, te, 0.0001);

		System.out.println("Linear-Gaussian TE validated");
	}
	
	public void testEmbedding() throws Exception {
		ArrayFileReader afr = new ArrayFileReader("demos/data/2coupledRandomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();

		int[] ks = {1,2,3,4};
		int[] tau_ks = {1,2,3,4};
		int[] ls = {1,2,3,4};
		int[] tau_ls = {1,2,3,4};
		int[] delays = {1,2,3,4};
		
		for (int kIndex = 0; kIndex < ks.length; kIndex++) {
			int k = ks[kIndex];
			for (int taukIndex = 0; taukIndex < tau_ks.length; taukIndex++) {
				int tau_k = tau_ks[taukIndex];
				for (int lIndex = 0; lIndex < ls.length; lIndex++) {
					int l = ls[lIndex];
					for (int taulIndex = 0; taulIndex < tau_ls.length; taulIndex++) {
						int tau_l = tau_ls[taulIndex];
						for (int delayIndex = 0; delayIndex < delays.length; delayIndex++) {
							int delay = delays[delayIndex];
							
							TransferEntropyCalculatorGaussian teCalc = new TransferEntropyCalculatorGaussian();
							teCalc.initialise(k, tau_k, l, tau_l, delay);
							teCalc.setObservations(MatrixUtils.selectColumn(data, 0), MatrixUtils.selectColumn(data, 1));
							
							// Compute the time index of the destination next step (indexed from 0)
							//  of the first complete observation tuple:
							int firstDestNextIndexIfDestPastLonger = (k-1)*tau_k + 1;
							int firstDestNextIndexIfSourcePastPlusDelayLonger = (l-1)*tau_l + 1 + (delay - 1);
							int firstDestNextIndex = Math.max(firstDestNextIndexIfDestPastLonger, firstDestNextIndexIfSourcePastPlusDelayLonger); 
							
							assertEquals(data.length -  firstDestNextIndex,
									teCalc.getNumObservations());
						}
					}
				}
			}
		}
	}
	
	public void testComputeStartEndPairsForSimpleEmbeddingTestData() throws Exception {
		
		boolean[] sourceValid = {true, true, true, false, true, false, true, true, true, true, true, true, true};
		boolean[] destValid =   {true, true, false, true, true, true, true, true, false, false, true, true, true};
		Vector<int[]> expected = new Vector<int[]>();
		expected.add(new int[] {0, 1});
		expected.add(new int[] {4, 5});
		expected.add(new int[] {6, 7});
		expected.add(new int[] {10, 12});
		
		int k = 1;
		TransferEntropyCalculatorGaussian teCalc =
				new TransferEntropyCalculatorGaussian();
		teCalc.initialise(k, 1, 1, 1, 1);
		Vector<int[]> timePairs = teCalc.computeStartAndEndTimePairs(sourceValid, destValid);
		// And compare to the existing 
		Vector<int[]> expectedTimePairs = computeStartAndEndTimePairs(sourceValid, destValid, k);
		
		assertEquals(expected.size(), timePairs.size());				
		assertEquals(expectedTimePairs.size(), timePairs.size());				
		for (int i = 0; i < expectedTimePairs.size(); i++) {
			int[] pair = timePairs.get(i);
			int[] expectedPairOldMethod = expectedTimePairs.get(i);
			int[] myExpectedPair = expected.get(i);
			assertEquals(myExpectedPair[0], pair[0]);
			assertEquals(myExpectedPair[1], pair[1]);
			assertEquals(expectedPairOldMethod[0], pair[0]);
			assertEquals(expectedPairOldMethod[1], pair[1]);
		}
	}

	public void testComputeStartEndPairsForEmbeddingTestDataK2() throws Exception {
		
		boolean[] sourceValid = {true, true, true, false, true, false, true, true, true, true, true, true, true, true};
		boolean[] destValid =   {true, true, false, true, true, true, true, true, false, false, true, true, true, true};
		Vector<int[]> expected = new Vector<int[]>();
		expected.add(new int[] {3, 5});
		expected.add(new int[] {5, 7});
		expected.add(new int[] {10, 13});
		
		int k = 2;
		TransferEntropyCalculatorGaussian teCalc =
				new TransferEntropyCalculatorGaussian();
		teCalc.initialise(k, 1, 1, 1, 1);
		Vector<int[]> timePairs = teCalc.computeStartAndEndTimePairs(sourceValid, destValid);
		
		assertEquals(expected.size(), timePairs.size());				
		for (int i = 0; i < timePairs.size(); i++) {
			int[] pair = timePairs.get(i);
			int[] myExpectedPair = expected.get(i);
			assertEquals(myExpectedPair[0], pair[0]);
			assertEquals(myExpectedPair[1], pair[1]);
		}
	}
	
	public void testComputeStartEndPairsForEmbeddingTestDataK2Delay() throws Exception {
		
		boolean[] sourceValid = {true, true, true, false, true, false, true, true, true, true, true, true, true, true};
		boolean[] destValid =   {true, true, false, true, true, true, true, true, false, false, true, true, true, true};
		Vector<int[]> expected = new Vector<int[]>();
		expected.add(new int[] {4, 6});
		expected.add(new int[] {10, 13});
		
		int k = 2;
		TransferEntropyCalculatorGaussian teCalc =
				new TransferEntropyCalculatorGaussian();
		teCalc.initialise(k, 1, 1, 1, 2);
		Vector<int[]> timePairs = teCalc.computeStartAndEndTimePairs(sourceValid, destValid);
		
		assertEquals(expected.size(), timePairs.size());				
		for (int i = 0; i < timePairs.size(); i++) {
			int[] pair = timePairs.get(i);
			int[] myExpectedPair = expected.get(i);
			assertEquals(myExpectedPair[0], pair[0]);
			assertEquals(myExpectedPair[1], pair[1]);
		}
	}

	public void testComputeStartEndPairsForEmbeddingTestDataK2LongerDelay() throws Exception {
		
		boolean[] sourceValid = {true, true, true, false, true, false, true, true, true, true, true, true, true, true};
		boolean[] destValid =   {true, true, false, true, true, true, true, true, false, false, true, true, true, true};
		Vector<int[]> expected = new Vector<int[]>();
		expected.add(new int[] {2, 5});
		expected.add(new int[] {4, 7});
		expected.add(new int[] {9, 13});
		
		int k = 2;
		TransferEntropyCalculatorGaussian teCalc =
				new TransferEntropyCalculatorGaussian();
		teCalc.initialise(k, 1, 1, 1, 3);
		Vector<int[]> timePairs = teCalc.computeStartAndEndTimePairs(sourceValid, destValid);
		
		assertEquals(expected.size(), timePairs.size());				
		for (int i = 0; i < timePairs.size(); i++) {
			int[] pair = timePairs.get(i);
			int[] myExpectedPair = expected.get(i);
			assertEquals(myExpectedPair[0], pair[0]);
			assertEquals(myExpectedPair[1], pair[1]);
		}
	}

	public void testComputeStartEndPairsForEmbeddingTestDataKLTausDelay() throws Exception {
		
		boolean[] sourceValid = {true, true, true, true, true, false, true, true, true, true, true, true, true, true};
		boolean[] destValid =   {true, true, false, true, true, true, true, true, true, false, true, true, true, true};
		
		Vector<int[]> expected = new Vector<int[]>();
		expected.add(new int[] {2, 6});
		expected.add(new int[] {9, 13});
		
		TransferEntropyCalculatorGaussian teCalc =
				new TransferEntropyCalculatorGaussian();
		teCalc.initialise(2, 2, 2, 2, 2);
		Vector<int[]> timePairs = teCalc.computeStartAndEndTimePairs(sourceValid, destValid);
		
		assertEquals(expected.size(), timePairs.size());				
		for (int i = 0; i < timePairs.size(); i++) {
			int[] pair = timePairs.get(i);
			int[] myExpectedPair = expected.get(i);
			assertEquals(myExpectedPair[0], pair[0]);
			assertEquals(myExpectedPair[1], pair[1]);
		}
	}

	public void testComputeStartEndPairsForSimpleEmbeddingRandomData() throws Exception {
		RandomGenerator rand = new RandomGenerator();
		
		// Generate random data on whether the sources and dests
		//  are valid, biasing towards valid
		int length = 1000;
		double[] randDoubles1 = rand.generateRandomData(length);
		double[] randDoubles2 = rand.generateRandomData(length);
		double pValid = 0.95;
		boolean[] sourceValid = new boolean[length];
		boolean[] destValid = new boolean[length];
		for (int t = 0; t < length; t++) {
			sourceValid[t] = randDoubles1[t] < pValid;
			destValid[t] = randDoubles2[t] < pValid;
		}
		
		int[] ks = {1,2,3,4};
		for (int kIndex = 0; kIndex < ks.length; kIndex++) {
			int k = ks[kIndex];
			TransferEntropyCalculatorGaussian teCalc =
					new TransferEntropyCalculatorGaussian();
			teCalc.initialise(k, 1, 1, 1, 1);
			Vector<int[]> timePairs = teCalc.computeStartAndEndTimePairs(sourceValid, destValid);
			// And compare to the existing 
			Vector<int[]> expectedTimePairs = computeStartAndEndTimePairs(sourceValid, destValid, k);
			
			assertEquals(expectedTimePairs.size(), timePairs.size());				
			for (int i = 0; i < expectedTimePairs.size(); i++) {
				int[] pair = timePairs.get(i);
				int[] expectedPair = expectedTimePairs.get(i);
				assertEquals(expectedPair[0], pair[0]);
				assertEquals(expectedPair[1], pair[1]);
			}
			System.out.printf("TransferEntropyViaCondMI.computeStartEndPairsForSimpleEmbedding correct for k=%d\n", k);
		}
	}

	/**
	 * This method copied from TransferEntropyCommon (since we may delete this soon)
	 *  to test simple cases of TransferEntropyCalculatorViaCondMutualInfo.computeStartAndEndTimePairs()
	 *  (i.e. where l, l_tau and delay are 1.)
	 * 
	 * @param sourceValid
	 * @param destValid
	 * @return
	 */
	public Vector<int[]> computeStartAndEndTimePairs(boolean[] sourceValid, boolean[] destValid, int k) {
		// Scan along the data avoiding invalid values
		int startTime = 0;
		int endTime = 0;
		boolean lookingForStart = true;
		Vector<int[]> startAndEndTimePairs = new Vector<int[]>();
		for (int t = 0; t < destValid.length; t++) {
			if (lookingForStart) {
				// Precondition: startTime holds a candidate start time
				if (destValid[t]) {
					// This point is OK at the destination
					if (t - startTime < k) {
						// We're still checking the past history only, so
						continue;
					} else {
						// We've got the full past history ok, so check the source also
						if (sourceValid[t - 1]) {
							// source is good to go also
							// set a candidate endTime
							endTime = t;
							lookingForStart = false;
							if (t == destValid.length - 1) {
								// we need to terminate now
								int[] timePair = new int[2];
								timePair[0] = startTime;
								timePair[1] = endTime;
								startAndEndTimePairs.add(timePair);
								// System.out.printf("t_s=%d, t_e=%d\n", startTime, endTime);
							}
						} else {
							// source was not good to go, so try moving along one time point
							startTime++;
						}
					}
				} else {
					// We need to keep looking.
					// Move the potential start time to the next point
					startTime = t + 1;
				}
			} else {
				// Precondition: startTime holds the start time for this set, 
				//  endTime holds a candidate end time
				// Check if we can include the current time step
				boolean terminateSequence = false;
				if (destValid[t] && sourceValid[t - 1]) {
					// We can extend
					endTime = t;
				} else {
					terminateSequence = true;
				}
				if (t == destValid.length - 1) {
					// we need to terminate the sequence anyway
					terminateSequence = true;
				}
				if (terminateSequence) {
					// This section is done
					int[] timePair = new int[2];
					timePair[0] = startTime;
					timePair[1] = endTime;
					startAndEndTimePairs.add(timePair);
					// System.out.printf("t_s=%d, t_e=%d\n", startTime, endTime);
					lookingForStart = true;
					if (!destValid[t]) {
						// The current destination observation broke our chain;
						//  so we need to start looking all over again:
						startTime = t + 1;
					} else {
						// The current source broke our chain (or we're at the
						//  end of the series anyway, so this doesn't matter);
						//  so we can keep the good destination history
						//  that we've built up here:
						startTime = t - k + 1;
					}
				}
			}
		}
		return startAndEndTimePairs;
	}
}
