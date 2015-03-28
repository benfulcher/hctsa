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

package infodynamics.measures.continuous;

import junit.framework.TestCase;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

public abstract class ConditionalTransferEntropyAbstractTester extends TestCase {

	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 * @param teCalc a pre-constructed ConditionalTransferEntropyCalculator object
	 * @param timeSteps number of time steps for the random data
	 * @param k dest history length for the TE calculator to use
	 */
	public void testLocalsAverageCorrectly(ConditionalTransferEntropyCalculator teCalc,
			int timeSteps, int k)
			throws Exception {
		
		teCalc.initialise(k, 1, 1);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[] sourceData = rg.generateNormalData(timeSteps,
				0, 1);
		double[] destData = rg.generateNormalData(timeSteps,
				0, 1);
		double[] condData = rg.generateNormalData(timeSteps,
				0, 1);
		
		teCalc.setObservations(sourceData, destData, condData);
		
		//teCalc.setDebug(true);
		double te = teCalc.computeAverageLocalOfObservations();
		//teCalc.setDebug(false);
		double[] teLocal = teCalc.computeLocalOfPreviousObservations();
		
		System.out.printf("Average was %.5f\n", te);

		assertEquals(te, MatrixUtils.mean(teLocal, k, timeSteps-k), 0.00001);
	}
	
	/**
	 * Confirm that significance testing doesn't alter the average that
	 * would be returned.
	 * 
	 * @param teCalc a pre-constructed ConditionalTransferEntropyCalculator object
	 * @param timeSteps number of time steps for the random data
	 * @param k history length for the TE calculator to use
	 * @throws Exception
	 */
	public void testComputeSignificanceDoesntAlterAverage(ConditionalTransferEntropyCalculator teCalc,
			int timeSteps, int k) throws Exception {
		
		teCalc.initialise(k, 1, 1);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[] sourceData = rg.generateNormalData(timeSteps,
				0, 1);
		double[] destData = rg.generateNormalData(timeSteps,
				0, 1);
		double[] condData = rg.generateNormalData(timeSteps,
				0, 1);
		
		teCalc.setObservations(sourceData, destData, condData);
		
		//teCalc.setDebug(true);
		double te = teCalc.computeAverageLocalOfObservations();
		//teCalc.setDebug(false);
		//double[] teLocal = teCalc.computeLocalOfPreviousObservations();
		
		System.out.printf("Average was %.5f\n", te);
		
		// Now look at statistical significance tests
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(
				timeSteps - k, 100);

		teCalc.computeSignificance(newOrderings);
		
		// And compute the average value again to check that it's consistent:
		for (int i = 0; i < 10; i++) {
			double averageCheck1 = teCalc.computeAverageLocalOfObservations();
			assertEquals(te, averageCheck1);
		}
	}

	/**
	 * Confirm that univariate method signature calls fail if the calculator
	 *  was not initialised for univariate conditional data.
	 * 
	 * @param teCalc
	 */
	public void testUnivariateCallFailsIfWrongInitialisation(ConditionalTransferEntropyCalculator teCalc) throws Exception {
		
		teCalc.initialise(1, 1, 1);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[] sourceData = rg.generateNormalData(10,
				0, 1);
		double[] destData = rg.generateNormalData(10,
				0, 1);
		double[] condData = rg.generateNormalData(10,
				0, 1);
		boolean gotException = false;
		try {
			teCalc.setObservations(sourceData, destData, condData);
		} catch (Exception e) {
			gotException = true;
		}
		assert(gotException);
	}
	
	/**
	 * Confirm the workings of the conditional TE calculator
	 *  by calculating pairwise TE with a long k, then computing
	 *  conditional TE with a shorter history length k, but
	 *  other conditional variables copying the values of those further past
	 *  values of the destination.
	 * 
	 * @param teCalc a pre-constructed TransferEntropyCalculator object
	 * @param condTeCalc a pre-constructed ConditionalTransferEntropyCalculator object
	 *   of the same estimator type
	 * @param timeSteps number of time steps for the random data
	 * @param k history length for the TE calculator to use
	 * @throws Exception
	 */
	public void testConditionalAgainstOrdinaryTE(
			TransferEntropyCalculator teCalc,
			ConditionalTransferEntropyCalculator condTeCalc,
			int timeSteps, int k) throws Exception {
		
		if (k < 2) {
			throw new Exception("Need k >= 2 for testConditionalAgainstOrdinaryTE");
		}
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[] sourceData = rg.generateNormalData(timeSteps,
				0, 1);
		double[] destData = rg.generateNormalData(timeSteps,
				0, 1);

		// First compute ordinary TE
		teCalc.initialise(k);
		teCalc.setObservations(sourceData, destData);
		double te = teCalc.computeAverageLocalOfObservations();
		System.out.printf("TE(k=%d): Average was %.5f\n", k, te);

		// Next compute conditional TE with some older
		//  parts of the destination as conditional variables
		//  instead of in the destination past.
		int[] condDims = {k-1};
		int[] condTaus = {1};
		int[] condDelays = {2};
		condTeCalc.initialise(1, 1, 1, 1, 1, condDims, condTaus, condDelays);
		// Don't need to extract the data ourselves:
		// double[][] condData =
		//		MatrixUtils.makeDelayEmbeddingVector(destData, k-1,
		//				k-2, destData.length-k+1); // Need an extra unused one here
		condTeCalc.setObservations(sourceData, destData, destData);
		double condTe = condTeCalc.computeAverageLocalOfObservations();
		System.out.printf("CondTE(k=%d): Average was %.5f\n", k, condTe);
		
		assertEquals(teCalc.getNumObservations(), condTeCalc.getNumObservations());
		assertEquals(te, condTe, 0.000000001);
	}
	
}
