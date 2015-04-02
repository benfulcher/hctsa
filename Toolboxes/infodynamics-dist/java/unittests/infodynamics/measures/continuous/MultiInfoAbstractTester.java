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

public abstract class MultiInfoAbstractTester extends TestCase {

	protected double lastResult = 0.0;
	
	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 * @param miCalc a pre-constructed MultiInfoCalculator object
	 * @param dimensions number of dimensions for the data to use
	 * @param timeSteps number of time steps for the random data
	 */
	public void testLocalsAverageCorrectly(MultiInfoCalculator miCalc,
			int dimensions, int timeSteps)
			throws Exception {
		
		miCalc.initialise(dimensions);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[][] data = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		
		miCalc.setObservations(data);
		
		//miCalc.setDebug(true);
		double mi = miCalc.computeAverageLocalOfObservations();
		lastResult = mi;
		//miCalc.setDebug(false);
		double[] miLocal = miCalc.computeLocalOfPreviousObservations();
		
		System.out.printf("Average was %.5f\n", mi);

		assertEquals(mi, MatrixUtils.mean(miLocal), 0.00001);
	}
	
	/**
	 * Confirm that significance testing doesn't alter the average that
	 * would be returned.
	 * 
	 * @param miCalc a pre-constructed MultiInfoCalculator object
	 * @param dimensions number of dimensions for the data to use
	 * @param timeSteps number of time steps for the random data
	 * @throws Exception
	 */
	public void testComputeSignificanceDoesntAlterAverage(MultiInfoCalculator miCalc,
			int dimensions, int timeSteps) throws Exception {
		
		miCalc.initialise(dimensions);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[][] data = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		
		miCalc.setObservations(data);
		
		//miCalc.setDebug(true);
		double mi = miCalc.computeAverageLocalOfObservations();
		//miCalc.setDebug(false);
		//double[] miLocal = miCalc.computeLocalOfPreviousObservations();
		
		System.out.printf("Average was %.5f\n", mi);
		
		// Now look at statistical significance tests
		EmpiricalMeasurementDistribution measDist = 
				miCalc.computeSignificance(2);
		// Make sure that (the first) surrogate TE does not
		//  match the actual TE (it could possibly match but with
		//  an incredibly low probability)
		assertFalse(mi == measDist.distribution[0]);
		
		// And compute the average value again to check that it's consistent:
		for (int i = 0; i < 10; i++) {
			double lastAverage = miCalc.getLastAverage();
			assertEquals(mi, lastAverage);
			double averageCheck1 = miCalc.computeAverageLocalOfObservations();
			assertEquals(mi, averageCheck1);
		}
	}

}
