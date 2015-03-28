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

import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;

public abstract class ConditionalMutualInfoMultiVariateAbstractTester
	extends TestCase {

	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 * @param condMiCalc a pre-constructed ConditionalMutualInfoCalculatorMultiVariate object
	 * @param dimensions number of dimensions for the source and dest data to use
	 * @param timeSteps number of time steps for the random data
	 */
	public void testLocalsAverageCorrectly(ConditionalMutualInfoCalculatorMultiVariate condMiCalc,
			int dimensions, int timeSteps)
			throws Exception {
		
		condMiCalc.initialise(dimensions, dimensions, dimensions);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[][] sourceData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] destData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] condData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		
		condMiCalc.setObservations(sourceData, destData, condData);
		
		//teCalc.setDebug(true);
		double condmi = condMiCalc.computeAverageLocalOfObservations();
		//miCalc.setDebug(false);
		double[] condMiLocal = condMiCalc.computeLocalOfPreviousObservations();
		
		System.out.printf("Average was %.5f\n", condmi);

		assertEquals(condmi, MatrixUtils.mean(condMiLocal), 0.00001);
	}
	
	/**
	 * Confirm that significance testing doesn't alter the average that
	 * would be returned.
	 * 
	 * @param condMiCalc a pre-constructed ConditionalMutualInfoCalculatorMultiVariate object
	 * @param dimensions number of dimensions for the source and dest data to use
	 * @param timeSteps number of time steps for the random data
	 * @throws Exception
	 */
	public void testComputeSignificanceDoesntAlterAverage(ConditionalMutualInfoCalculatorMultiVariate condMiCalc,
			int dimensions, int timeSteps) throws Exception {
		
		condMiCalc.initialise(dimensions, dimensions, dimensions);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[][] sourceData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] destData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] condData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		
		condMiCalc.setObservations(sourceData, destData, condData);
		
		//condMiCalc.setDebug(true);
		double condMi = condMiCalc.computeAverageLocalOfObservations();
		//condMiCalc.setDebug(false);
		//double[] condMiLocal = miCalc.computeLocalOfPreviousObservations();
		
		System.out.printf("Average was %.5f\n", condMi);
		
		// Now look at statistical significance tests
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(
				timeSteps, 100);
		
		// Compute significance for permuting first variable 
		EmpiricalMeasurementDistribution measDist =
				condMiCalc.computeSignificance(1, newOrderings);
		// The actual MI should be different to a surrogate (it's possible
		//  but exceedingly unlikely that they would be equal).
		assertFalse(condMi == measDist.distribution[0]);
		
		// And compute the average value again to check that it's consistent:
		for (int i = 0; i < 10; i++) {
			double lastAverage = condMiCalc.getLastAverage();
			assertEquals(condMi, lastAverage);
			double averageCheck1 = condMiCalc.computeAverageLocalOfObservations();
			assertEquals(condMi, averageCheck1);
		}

		// Compute significance for permuting second variable 
		condMiCalc.computeSignificance(2, newOrderings);
		
		// And compute the average value again to check that it's consistent:
		for (int i = 0; i < 10; i++) {
			double lastAverage = condMiCalc.getLastAverage();
			assertEquals(condMi, lastAverage);
			double averageCheck1 = condMiCalc.computeAverageLocalOfObservations();
			assertEquals(condMi, averageCheck1);
		}
	}

}
