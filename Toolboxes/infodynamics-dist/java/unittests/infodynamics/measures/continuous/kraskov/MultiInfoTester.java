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

package infodynamics.measures.continuous.kraskov;

import infodynamics.measures.continuous.MultiInfoAbstractTester;
import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;

public class MultiInfoTester extends MultiInfoAbstractTester {

	protected String NUM_THREADS_TO_USE_DEFAULT = MultiInfoCalculatorKraskov.USE_ALL_THREADS;
	protected String NUM_THREADS_TO_USE = NUM_THREADS_TO_USE_DEFAULT;
	
	/**
	 * Utility function to create a calculator for the given algorithm number
	 * 
	 * @param algNumber
	 * @return
	 */
	public MultiInfoCalculatorKraskov getNewCalc(int algNumber) {
		MultiInfoCalculatorKraskov miCalc = null;
		if (algNumber == 1) {
			miCalc = new MultiInfoCalculatorKraskov1();
		} else if (algNumber == 2) {
			miCalc = new MultiInfoCalculatorKraskov2();
		}
		return miCalc;
	}
	
	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 * 
	 */
	public void checkLocalsAverageCorrectly(int algNumber, String numThreads) throws Exception {
		
		MultiInfoCalculatorKraskov miCalc = getNewCalc(algNumber);
		
		String kraskov_K = "4";
		
		miCalc.setProperty(
				MultiInfoCalculatorKraskov.PROP_K,
				kraskov_K);
		miCalc.setProperty(
				MultiInfoCalculatorKraskov.PROP_NUM_THREADS,
				numThreads);

		super.testLocalsAverageCorrectly(miCalc, 2, 10000);
	}
	public void testLocalsAverageCorrectly() throws Exception {
		checkLocalsAverageCorrectly(1, NUM_THREADS_TO_USE);
		checkLocalsAverageCorrectly(2, NUM_THREADS_TO_USE);
	}
	
	/**
	 * Confirm that significance testing doesn't alter the average that
	 * would be returned.
	 * 
	 * @throws Exception
	 */
	public void checkComputeSignificanceDoesntAlterAverage(int algNumber) throws Exception {
		
		MultiInfoCalculatorKraskov miCalc = getNewCalc(algNumber);
		
		String kraskov_K = "4";
		
		miCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);
		miCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
				NUM_THREADS_TO_USE);

		super.testComputeSignificanceDoesntAlterAverage(miCalc, 2, 100);
	}
	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		checkComputeSignificanceDoesntAlterAverage(1);
		checkComputeSignificanceDoesntAlterAverage(2);
	}
	
	/**
	 * Utility function to run Kraskov MI for data with known results
	 * 
	 * @param var1
	 * @param kNNs array of Kraskov k nearest neighbours parameter to check
	 * @param expectedResults array of expected results for each k
	 */
	protected void checkMIForGivenData(double[][] data, 
			int[] kNNs, double[] expectedResults) throws Exception {
				
		// The Kraskov MILCA toolkit MIhigherdim executable 
		//  uses algorithm 2 by default (this is what it means by rectangular):
		MultiInfoCalculatorKraskov miCalc = getNewCalc(2);
		
		for (int kIndex = 0; kIndex < kNNs.length; kIndex++) {
			int k = kNNs[kIndex];
			miCalc.setProperty(
					MutualInfoCalculatorMultiVariateKraskov.PROP_K,
					Integer.toString(k));
			miCalc.setProperty(
					MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
					NUM_THREADS_TO_USE);
			// No longer need to set this property as it's set by default:
			//miCalc.setProperty(MutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE,
			//		EuclideanUtils.NORM_MAX_NORM_STRING);
			miCalc.setProperty(MultiInfoCalculatorKraskov.PROP_ADD_NOISE, "0"); // Need consistency of results for unit test
			miCalc.initialise(data[0].length);
			miCalc.setObservations(data);
			miCalc.setDebug(true);
			double mi = miCalc.computeAverageLocalOfObservations();
			miCalc.setDebug(false);
			
			System.out.printf("k=%d: Average Multi-info %.8f (expected %.8f)\n",
					k, mi, expectedResults[kIndex]);
			// Dropping required accuracy by one order of magnitude, due
			//  to faster but slightly less accurate digamma estimator change
			// 0.0000001 is fine for all but last test, so dropping again
			//  to 0.000001
			assertEquals(expectedResults[kIndex], mi, 0.000001);			
		}
	}
	
	/**
	 * Test the computed Multi info for 2 variables (i.e. should be a regular MI!)
	 *  against that calculated by Kraskov's own MILCA
	 * tool on the same data.
	 * 
	 * To run Kraskov's tool (http://www.klab.caltech.edu/~kraskov/MILCA/) for this 
	 * data, run:
	 * ./MIxnyn <dataFile> 1 1 3000 <kNearestNeighbours> 0
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testUnivariateMIforRandomVariablesFromFile() throws Exception {
		
		// Test set 1:
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/2randomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {1, 2, 3, 4, 5, 6, 10, 15};
		// Expected values from Kraskov's MILCA toolkit:
		double[] expectedFromMILCA = {-0.05294175, -0.03944338, -0.02190217,
				0.00120807, -0.00924771, -0.00316402, -0.00778205, -0.00565778};
		
		System.out.println("Kraskov comparison 1 - univariate random data 1");
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0, 1}),
				kNNs, expectedFromMILCA);
		
		//------------------
		// Test set 2:
		
		// We'll just take the first two columns from this data set
		afr = new ArrayFileReader("demos/data/4randomCols-1.txt");
		data = afr.getDouble2DMatrix();
		
		// Expected values from Kraskov's MILCA toolkit:
		double[] expectedFromMILCA_2 = {-0.04614525, -0.00861460, -0.00164540,
				-0.01130354, -0.01339670, -0.00964035, -0.00237072, -0.00096891};
		
		System.out.println("Kraskov comparison 2 - univariate random data 2");
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0, 1}),
				kNNs, expectedFromMILCA_2);

	}

	/**
	 * Test the computed multivariate multi-info against that calculated by Kraskov's own MILCA
	 * tool on the same data.
	 * 
	 * To run Kraskov's tool (http://www.klab.caltech.edu/~kraskov/MILCA/) for this 
	 * data, run:
	 * ./MIhigherdim <dataFile> 4 1 1 3000 <kNearestNeighbours> 0
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testMultivariateMIforRandomVariablesFromFile() throws Exception {
		
		// Test set 3:
		
		// We'll just take the first two columns from this data set
		ArrayFileReader afr = new ArrayFileReader("demos/data/4randomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {1, 2, 3, 4, 5, 6, 10, 15};
		// Expected values from Kraskov's MILCA toolkit:
		double[] expectedFromMILCA_2 = {0.03229833, -0.01146200, -0.00691358,
				0.00002149, -0.01056322, -0.01482730, -0.01223885, -0.01461794};
		
		System.out.println("Kraskov comparison 3 - multivariate random data 1");
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0, 1, 2, 3}),
				kNNs, expectedFromMILCA_2);

	}

	/**
	 * Test the computed multivariate multi-info against that calculated by Kraskov's own MILCA
	 * tool on the same data.
	 * 
	 * To run Kraskov's tool (http://www.klab.caltech.edu/~kraskov/MILCA/) for this 
	 * data, run:
	 * ./MIhigherdim <dataFile> 4 1 1 3000 <kNearestNeighbours> 0
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testMultivariateMIVariousNumThreads() throws Exception {
		
		// Test set 3:
		
		// We'll just take the first two columns from this data set
		ArrayFileReader afr = new ArrayFileReader("demos/data/4randomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {3, 4};
		// Expected values from Kraskov's MILCA toolkit:
		double[] expectedFromMILCA_2 = {-0.00691358,
				0.00002149};
		
		System.out.println("Kraskov comparison 3a - single threaded");
		NUM_THREADS_TO_USE = "1";
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0, 1, 2, 3}),
				kNNs, expectedFromMILCA_2);
		System.out.println("Kraskov comparison 3b - dual threaded");
		NUM_THREADS_TO_USE = "2";
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0, 1, 2, 3}),
				kNNs, expectedFromMILCA_2);
		NUM_THREADS_TO_USE = NUM_THREADS_TO_USE_DEFAULT;
	}

	/**
	 * Test the computed multivariate MI against that calculated by Kraskov's own MILCA
	 * tool on the same data.
	 * 
	 * To run Kraskov's tool (http://www.klab.caltech.edu/~kraskov/MILCA/) for this 
	 * data, run:
	 * ./MIhigherdim <dataFile> 4 1 1 3000 <kNearestNeighbours> 0
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testMultivariateMIforDependentVariablesFromFile() throws Exception {
		
		// Test set 6:
		
		// We'll just take the first two columns from this data set
		ArrayFileReader afr = new ArrayFileReader("demos/data/4ColsPairedDirectDependence-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {1, 2, 3, 4, 5, 6, 10, 15};
		// Expected values from Kraskov's MILCA toolkit:
		double[] expectedFromMILCA_2 = {8.44056282, 7.69813699, 7.26909347, 
				6.97095249, 6.73728113, 6.53105867, 5.96391264, 5.51627278};
		
		System.out.println("Kraskov comparison 6 - multivariate dependent data 1");
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0, 1, 2, 3}),
				kNNs, expectedFromMILCA_2);

	}

	/**
	 * Test the computed multivariate MI against that calculated by Kraskov's own MILCA
	 * tool on the same data.
	 * 
	 * To run Kraskov's tool (http://www.klab.caltech.edu/~kraskov/MILCA/) for this 
	 * data, run:
	 * ./MIhigherdim <dataFile> 4 1 1 3000 <kNearestNeighbours> 0
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testMultivariateMIforNoisyDependentVariablesFromFile() throws Exception {
		
		// Test set 7:
		
		// We'll just take the first two columns from this data set
		ArrayFileReader afr = new ArrayFileReader("demos/data/4ColsPairedNoisyDependence-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {1, 2, 3, 4, 5, 6, 10, 15};
		// Expected values from Kraskov's MILCA toolkit:
		double[] expectedFromMILCA_2 = {0.31900665, 0.37304998, 0.37213228, 
				0.37982388, 0.37304217, 0.36802502, 0.36353436, 0.35095074};
		
		System.out.println("Kraskov comparison 7 - multivariate dependent data 1");
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0, 1, 2, 3}),
				kNNs, expectedFromMILCA_2);

	}

	/**
	 * Test the computed multivariate MI against that calculated by Kraskov's own MILCA
	 * tool on the same data.
	 * 
	 * To run Kraskov's tool (http://www.klab.caltech.edu/~kraskov/MILCA/) for this 
	 * data, run:
	 * ./MIhigherdim <dataFile> 10 1 1 10000 <kNearestNeighbours> 0
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testMultivariateMIforRandomGaussianVariablesFromFile() throws Exception {
		
		// Test set 8:
		
		// We'll take the columns from this data set
		ArrayFileReader afr = new ArrayFileReader("demos/data/10ColsRandomGaussian-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {1, 2, 4, 10, 15};
		// Expected values from Kraskov's MILCA toolkit:
		double[] expectedFromMILCA_2 = {0.00932984, 0.00662195, 0.01697033,
				0.00397984, 0.00212609};
		
		System.out.println("Kraskov comparison 8 - multivariate uncorrelated Gaussian data 1");
		checkMIForGivenData(MatrixUtils.selectColumns(data,
								new int[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}),
				kNNs, expectedFromMILCA_2);

	}
}
