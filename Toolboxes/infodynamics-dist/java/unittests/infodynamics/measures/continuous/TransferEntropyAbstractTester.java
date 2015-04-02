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

public abstract class TransferEntropyAbstractTester extends TestCase {

	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 * @param teCalc a pre-constructed TransferEntropyCalculator object
	 * @param timeSteps number of time steps for the random data
	 * @param k dest history length for the TE calculator to use
	 */
	public void testLocalsAverageCorrectly(TransferEntropyCalculator teCalc,
			int timeSteps, int k)
			throws Exception {
		
		System.out.println("Testing locals average correctly");
		
		teCalc.initialise(k);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[] sourceData = rg.generateNormalData(timeSteps,
				0, 1);
		double[] destData = rg.generateNormalData(timeSteps,
				0, 1);
		
		teCalc.setObservations(sourceData, destData);
		
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
	 * @param teCalc a pre-constructed TransferEntropyCalculator object
	 * @param timeSteps number of time steps for the random data
	 * @param k history length for the TE calculator to use
	 * @throws Exception
	 */
	public void testComputeSignificanceDoesntAlterAverage(TransferEntropyCalculator teCalc,
			int timeSteps, int k) throws Exception {
		
		System.out.printf("Testing significance doesn't alter average");
		
		teCalc.initialise(k);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[] sourceData = rg.generateNormalData(timeSteps,
				0, 1);
		double[] destData = rg.generateNormalData(timeSteps,
				0, 1);
		
		teCalc.setObservations(sourceData, destData);
		
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

}
