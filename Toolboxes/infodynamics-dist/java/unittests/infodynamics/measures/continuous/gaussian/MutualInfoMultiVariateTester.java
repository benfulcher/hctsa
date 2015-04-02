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

import junit.framework.TestCase;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.MatrixUtilsTest;
import infodynamics.utils.RandomGenerator;

public class MutualInfoMultiVariateTester extends TestCase {

	public void testCovarianceDoesntMatchDimensions() throws Exception {
		MutualInfoCalculatorMultiVariateGaussian miCalc =
				new MutualInfoCalculatorMultiVariateGaussian();
		miCalc.initialise(1, 1);
		boolean caughtException = false;
		// Check that we catch a covariance matrix which doesn't have
		//  the right number of rows
		try {
			miCalc.setCovariance(new double[][] {{2,1,0.5}, {1,2,0.5}, {0.5,0.5,2}}, 1);
		} catch (Exception e) {
			caughtException = true;
		}
		assertTrue(caughtException);
		// Check that we catch a covariance matrix which isn't square
		caughtException = false;
		try {
			miCalc.setCovariance(new double[][] {{2,1}, {1,2,0.5}}, 1);
		} catch (Exception e) {
			caughtException = true;
		}
		assertTrue(caughtException);
		// Check that we catch a covariance matrix which isn't symmetric
		caughtException = false;
		try {
			miCalc.setCovariance(new double[][] {{2,1}, {1.000001,2}}, 1);
		} catch (Exception e) {
			caughtException = true;
		}
		assertTrue(caughtException);
		// Check that no exception is thrown for an ok covariance matrix
		caughtException = false;
		double[][] goodCovariance = new double[][] {{2,1}, {1,2}};
		try {
			miCalc.setCovariance(goodCovariance, 93);
		} catch (Exception e) {
			caughtException = true;
		}
		assertFalse(caughtException);
		// and that this covariance has been set (by verifying the Cholesky
		// decomposition of it is stored):
		MatrixUtilsTest.checkMatrix(MatrixUtils.CholeskyDecomposition(goodCovariance),
				miCalc.L, 0.00001);
	}
	
	public void testAnalyticMatchesCouplingValue() throws Exception {
		double[][] covarianceMatrix = {{1.0, 0}, {0, 1.0}};
		
		MutualInfoCalculatorMultiVariateGaussian miCalc =
				new MutualInfoCalculatorMultiVariateGaussian();
		
		for (double covar = 0.0; covar < 1.0; covar += 0.01) {
			// Insert the new covariance into the matrix:
			covarianceMatrix[0][1] = covar;
			covarianceMatrix[1][0] = covar;

			miCalc.initialise(1, 1);
			miCalc.setCovariance(covarianceMatrix, 1);
			assertEquals(-0.5 * Math.log(1.0 - covar*covar),
					miCalc.computeAverageLocalOfObservations(), 0.00000000001);
		}
	}
	
	public void testMIfromSuppliedCovariance() throws Exception {
		MutualInfoCalculatorMultiVariateGaussian miCalc =
				new MutualInfoCalculatorMultiVariateGaussian();

		double[][] covarianceMatrix = {{5, 3}, {3, 4}}; // det is 11
		miCalc.initialise(1, 1);
		miCalc.setCovariance(covarianceMatrix, 100);
		assertEquals(0.5 * Math.log(20.0 / 11.0),
				miCalc.computeAverageLocalOfObservations(), 0.00000000001);
		
		double[][] covarianceMatrix2 = {{5, 3, 1}, {3, 4, 1.5}, {1, 1.5, 2}}; // det is 15.75
		miCalc.initialise(2, 1);
		miCalc.setCovariance(covarianceMatrix2, 100);
		assertEquals(0.5 * Math.log(11.0 * 2.0 / 15.75), // marginal dets are 11 and 2
				miCalc.computeAverageLocalOfObservations(), 0.00000000001);

	}
	
	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		
		MutualInfoCalculatorMultiVariateGaussian miCalc =
				new MutualInfoCalculatorMultiVariateGaussian();
		
		int dimensions = 2;
		int timeSteps = 100;
		
		miCalc.initialise(dimensions, dimensions);
		
		// generate some random data
		RandomGenerator rg = new RandomGenerator();
		double[][] sourceData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		double[][] destData = rg.generateNormalData(timeSteps, dimensions,
				0, 1);
		
		miCalc.setObservations(sourceData, destData);
		
		//miCalc.setDebug(true);
		double mi = miCalc.computeAverageLocalOfObservations();
		//miCalc.setDebug(false);
		//double[] miLocal = miCalc.computeLocalOfPreviousObservations();
		
		System.out.printf("Average was %.5f\n", mi);
		
		// Now look at statistical significance tests
		int[][] newOrderings = rg.generateDistinctRandomPerturbations(
				timeSteps, 100);

		EmpiricalMeasurementDistribution measDist =
				miCalc.computeSignificance(newOrderings);
		
		System.out.printf("pValue of sig test was %.3f\n", measDist.pValue);
		
		// And compute the average value again to check that it's consistent:
		for (int i = 0; i < 10; i++) {
			double averageCheck1 = miCalc.computeAverageLocalOfObservations();
			assertEquals(mi, averageCheck1);
		}
	}
}
