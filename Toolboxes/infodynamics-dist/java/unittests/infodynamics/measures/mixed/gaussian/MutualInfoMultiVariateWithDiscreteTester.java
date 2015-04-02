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

package infodynamics.measures.mixed.gaussian;

import junit.framework.TestCase;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

public class MutualInfoMultiVariateWithDiscreteTester extends TestCase {

	int DATA_LENGTH = 100;
	int TEST_DIMENSIONS = 2;
	int TEST_BASE = 2;

	public void testSetObservationsDataIntegrity() throws Exception {	
		
		MutualInfoCalculatorMultiVariateWithDiscreteGaussian miCalc =
				new MutualInfoCalculatorMultiVariateWithDiscreteGaussian();
		miCalc.initialise(TEST_DIMENSIONS, TEST_BASE);
		
		// Generate proper and incorrect data sets
		RandomGenerator rg = new RandomGenerator();
		double[][] contDataTooManyVariables = rg.generateNormalData(DATA_LENGTH,
				TEST_DIMENSIONS + 1, 0, 1);
		double[][] contDataCorrect = rg.generateNormalData(DATA_LENGTH, TEST_DIMENSIONS,
				0, 1);
		int[] discDataCorrect = rg.generateRandomInts(DATA_LENGTH, TEST_BASE);
		int[] discDataValuesOutsideRange = rg.generateRandomInts(100, TEST_BASE + 1);
		int[] discDataWrongLength = rg.generateRandomInts(DATA_LENGTH - 1, TEST_BASE);

		// Check that we catch continuous data which doesn't have
		//  the right number of variables
		boolean caughtException = false;
		try {
			miCalc.setObservations(contDataTooManyVariables, discDataCorrect);
		} catch (Exception e) {
			caughtException = true;
		}
		assertTrue(caughtException);
		
		// Check that we catch discrete data outside the range:
		caughtException = false;
		try {
			miCalc.setObservations(contDataCorrect, discDataValuesOutsideRange);
		} catch (Exception e) {
			caughtException = true;
		}
		assertTrue(caughtException);
		
		// Check that we catch mismatched data lengths:
		caughtException = false;
		try {
			miCalc.setObservations(contDataCorrect, discDataWrongLength);
		} catch (Exception e) {
			caughtException = true;
		}
		assertTrue(caughtException);

		// No observations should have been set yet
		assertEquals(0, miCalc.getNumObservations());
		
		// Check that no exception is thrown for ok data:
		caughtException = false;
		try {
			miCalc.setObservations(contDataCorrect, discDataCorrect);
		} catch (Exception e) {
			caughtException = true;
			e.printStackTrace();
		}
		assertFalse(caughtException);
		
		// and that the observations have now been set:
		assertEquals(DATA_LENGTH, miCalc.getNumObservations());
	}
		
	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		
		MutualInfoCalculatorMultiVariateWithDiscreteGaussian miCalc =
				new MutualInfoCalculatorMultiVariateWithDiscreteGaussian();
		miCalc.initialise(TEST_DIMENSIONS, TEST_BASE);
		
		// Generate proper data sets
		RandomGenerator rg = new RandomGenerator();
		double[][] contDataCorrect = rg.generateNormalData(DATA_LENGTH, TEST_DIMENSIONS,
				0, 1);
		int[] discDataCorrect = rg.generateRandomInts(DATA_LENGTH, TEST_BASE);

		miCalc.setObservations(contDataCorrect, discDataCorrect);

		//miCalc.setDebug(true);
		double mi = miCalc.computeAverageLocalOfObservations();
		//miCalc.setDebug(false);
		//double[] miLocal = miCalc.computeLocalOfPreviousObservations();
		
		System.out.printf("Average was %.5f\n", mi);
		
		EmpiricalMeasurementDistribution measDist =
				miCalc.computeSignificance(100);
		
		System.out.printf("pValue of sig test was %.3f\n", measDist.pValue);
		
		// And compute the average value again to check that it's consistent:
		for (int i = 0; i < 10; i++) {
			double averageCheck1 = miCalc.computeAverageLocalOfObservations();
			assertEquals(mi, averageCheck1);
		}
	}
	
	public void testObservationsRequiredBeforeStatTest() throws Exception {
		MutualInfoCalculatorMultiVariateWithDiscreteGaussian miCalc =
				new MutualInfoCalculatorMultiVariateWithDiscreteGaussian();
		miCalc.initialise(TEST_DIMENSIONS, TEST_BASE);

		boolean caughtException = false;
		try {
			miCalc.computeSignificance(100);
		} catch (Exception e) {
			caughtException = true;
		}
		assertTrue(caughtException);
	}
}
