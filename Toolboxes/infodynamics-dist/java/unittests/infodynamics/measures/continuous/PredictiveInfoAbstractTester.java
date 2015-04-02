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

public abstract class PredictiveInfoAbstractTester extends TestCase {

	protected double lastResult = 0.0;
	
	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 * @param piCalc a pre-constructed PredictiveInfoCalculator object
	 * @param k embedding length for past and future to use
	 * @param timeSteps number of time steps for the random data
	 */
	public void testLocalsAverageCorrectly(PredictiveInfoCalculator piCalc,
			int k, int timeSteps)
			throws Exception {
		
		piCalc.initialise(k);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[] data = rg.generateNormalData(timeSteps, 
				0, 1);
		
		piCalc.setObservations(data);
		
		//piCalc.setDebug(true);
		double pi = piCalc.computeAverageLocalOfObservations();
		lastResult = pi;
		//piCalc.setDebug(false);
		double[] piLocal = piCalc.computeLocalOfPreviousObservations();
		
		System.out.printf("Average was %.5f\n", pi);

		assertEquals(pi, MatrixUtils.mean(piLocal), 0.00001);
	}
	
	/**
	 * Confirm that significance testing doesn't alter the average that
	 * would be returned.
	 * 
	 * @param piCalc a pre-constructed PredictiveInfoCalculator object
	 * @param k embedding length for past and future to use
	 * @param timeSteps number of time steps for the random data
	 * @throws Exception
	 */
	public void testComputeSignificanceDoesntAlterAverage(PredictiveInfoCalculator piCalc,
			int k, int timeSteps) throws Exception {
		
		piCalc.initialise(k);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[] data = rg.generateNormalData(timeSteps,
				0, 1);
		
		piCalc.setObservations(data);
		
		//piCalc.setDebug(true);
		double pi = piCalc.computeAverageLocalOfObservations();
		//piCalc.setDebug(false);
		//double[] piLocal = piCalc.computeLocalOfPreviousObservations();
		
		System.out.printf("Average was %.5f\n", pi);
		
		// Now look at statistical significance tests
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(
				timeSteps-(2*k-1), 2);

		EmpiricalMeasurementDistribution measDist = 
				piCalc.computeSignificance(newOrderings);
		// Make sure that (the first) surrogate TE does not
		//  match the actual TE (it could possibly match but with
		//  an incredibly low probability)
		assertFalse(pi == measDist.distribution[0]);
		
		// And compute the average value again to check that it's consistent:
		for (int i = 0; i < 10; i++) {
			double lastAverage = piCalc.getLastAverage();
			assertEquals(pi, lastAverage);
			double averageCheck1 = piCalc.computeAverageLocalOfObservations();
			assertEquals(pi, averageCheck1);
		}
	}

}
