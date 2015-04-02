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

import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;

public class MutualInfoMultiVariateTester
	extends infodynamics.measures.continuous.MutualInfoMultiVariateAbstractTester {

	protected String NUM_THREADS_TO_USE_DEFAULT = MutualInfoCalculatorMultiVariateKraskov.USE_ALL_THREADS;
	protected String NUM_THREADS_TO_USE = NUM_THREADS_TO_USE_DEFAULT;
	
	/**
	 * Utility function to create a calculator for the given algorithm number
	 * 
	 * @param algNumber
	 * @return
	 */
	public MutualInfoCalculatorMultiVariateKraskov getNewCalc(int algNumber) {
		MutualInfoCalculatorMultiVariateKraskov miCalc = null;
		if (algNumber == 1) {
			miCalc = new MutualInfoCalculatorMultiVariateKraskov1();
		} else if (algNumber == 2) {
			miCalc = new MutualInfoCalculatorMultiVariateKraskov2();
		}
		return miCalc;
	}
	
	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 * 
	 */
	public void checkLocalsAverageCorrectly(int algNumber, String numThreads) throws Exception {
		
		MutualInfoCalculatorMultiVariateKraskov miCalc = getNewCalc(algNumber);
		
		String kraskov_K = "4";
		
		miCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);
		miCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
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
		
		MutualInfoCalculatorMultiVariateKraskov miCalc = getNewCalc(algNumber);
		
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
	 * @param var2
	 * @param kNNs array of Kraskov k nearest neighbours parameter to check
	 * @param expectedResults array of expected results for each k
	 */
	protected void checkMIForGivenData(double[][] var1, double[][] var2,
			int[] kNNs, double[] expectedResults) throws Exception {
				
		// The Kraskov MILCA toolkit MIhigherdim executable 
		//  uses algorithm 2 by default (this is what it means by rectangular):
		MutualInfoCalculatorMultiVariateKraskov miCalc = getNewCalc(2);
		
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
			miCalc.initialise(var1[0].length, var2[0].length);
			miCalc.setObservations(var1, var2);
			miCalc.setDebug(true);
			double mi = miCalc.computeAverageLocalOfObservations();
			miCalc.setDebug(false);
			
			System.out.printf("k=%d: Average MI %.8f (expected %.8f)\n",
					k, mi, expectedResults[kIndex]);
			// Dropping required accuracy by one order of magnitude, due
			//  to faster but slightly less accurate digamma estimator change
			assertEquals(expectedResults[kIndex], mi, 0.0000001);			
		}
	}
	
	/**
	 * Test the computed univariate MI against that calculated by Kraskov's own MILCA
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
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
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
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				kNNs, expectedFromMILCA_2);

	}

	/**
	 * Test the computed multivariate MI against that calculated by Kraskov's own MILCA
	 * tool on the same data.
	 * 
	 * To run Kraskov's tool (http://www.klab.caltech.edu/~kraskov/MILCA/) for this 
	 * data, run:
	 * ./MIxnyn <dataFile> 2 2 3000 <kNearestNeighbours> 0
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
		double[] expectedFromMILCA_2 = {0.02886644, 0.01071634, 0.00186857,
				-0.00377259, -0.00634851, -0.00863725, -0.01058087, -0.01106348};
		
		System.out.println("Kraskov comparison 3 - multivariate random data 1");
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0, 1}),
				MatrixUtils.selectColumns(data, new int[] {2, 3}),
				kNNs, expectedFromMILCA_2);

	}

	/**
	 * Test the computed multivariate MI against that calculated by Kraskov's own MILCA
	 * tool on the same data, using various numbers of threads
	 * 
	 * To run Kraskov's tool (http://www.klab.caltech.edu/~kraskov/MILCA/) for this 
	 * data, run:
	 * ./MIxnyn <dataFile> 2 2 3000 <kNearestNeighbours> 0
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testMultivariateMIVariousNumThreads() throws Exception {
		
		// Test set 3a and 3b:
		
		// We'll just take the first two columns from this data set
		ArrayFileReader afr = new ArrayFileReader("demos/data/4randomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {3, 4};
		// Expected values from Kraskov's MILCA toolkit:
		double[] expectedFromMILCA_2 = {0.00186857,
				-0.00377259};
		
		System.out.println("Kraskov comparison 3a - single threaded");
		NUM_THREADS_TO_USE = "1";
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0, 1}),
				MatrixUtils.selectColumns(data, new int[] {2, 3}),
				kNNs, expectedFromMILCA_2);
		System.out.println("Kraskov comparison 3b - dual threaded");
		NUM_THREADS_TO_USE = "2";
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0, 1}),
				MatrixUtils.selectColumns(data, new int[] {2, 3}),
				kNNs, expectedFromMILCA_2);
		NUM_THREADS_TO_USE = NUM_THREADS_TO_USE_DEFAULT;

	}

	/**
	 * Test the computed multivariate MI against that calculated by Kraskov's own MILCA
	 * tool on the same data.
	 * 
	 * To run Kraskov's tool (http://www.klab.caltech.edu/~kraskov/MILCA/) for this 
	 * data, run:
	 * ./MIxnyn <dataFile> 1 3 3000 <kNearestNeighbours> 0
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testImbalancedMultivariateMIforRandomVariablesFromFile() throws Exception {
		
		// Test set 4:
		
		// We'll take MI from first column to the next 3:
		ArrayFileReader afr = new ArrayFileReader("demos/data/4randomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {1, 2, 3, 4, 5, 6, 10, 15};
		// Expected values from Kraskov's MILCA toolkit:
		double[] expectedFromMILCA = {0.02473475, 0.00404451, -0.00454679,
				-0.00737512, -0.00464896, -0.00610772, -0.00881741, -0.01306668};
		
		System.out.println("Kraskov comparison 4 - multivariate random data 2 (1 var to 3 vars)");
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1, 2, 3}),
				kNNs, expectedFromMILCA);
	}

	/**
	 * Test the computed multivariate MI against that calculated by Kraskov's own MILCA
	 * tool on the same data.
	 * This also tests for multithreading with residuals, assuming
	 *  we're running on a 4 processor machine
	 * 
	 * To run Kraskov's tool (http://www.klab.caltech.edu/~kraskov/MILCA/) for this 
	 * data, run:
	 * ./MIxnyn <dataFile> 2 2 3030 <kNearestNeighbours> 0
	 * where the file has the first 30 rows repeated.
	 * 
	 * Kraskov et al recommend that a small amount of noise should be 
	 * added to the data to avoid issues with repeated scores; this
	 * can be done in our toolkit by setting the relevant property
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testMultivariateMIforRandomVariablesRepeatedDataFromFile() throws Exception {
		
		// Test set 5:
		
		// We'll just take the first two columns from this data set
		ArrayFileReader afr = new ArrayFileReader("demos/data/4randomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		double[][] data2 = new double[data.length + 30][data[0].length];
		for (int r = 0; r < data.length; r++) {
			for (int c = 0; c < data[r].length; c++) {
				data2[r][c] = data[r][c];
			}
		}
		// Repeat the first 30 rows:
		for (int r = 0; r < 30; r++) {
			for (int c = 0; c < data[r].length; c++) {
				data2[r+data.length][c] = data[r][c];
			}			
		}
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {1, 2, 3, 4, 5, 6, 10, 15};
		// Expected values from Kraskov's MILCA toolkit:
		double[] expectedFromMILCA_2 = {0.16846374, 0.04091779, 0.02069109,
				0.00700680, 0.00121768, -0.00134164, -0.00870685, -0.00966508};
		
		System.out.println("Kraskov comparison 5 - multivariate random data 1 with 30 repeated rows");
		checkMIForGivenData(MatrixUtils.selectColumns(data2, new int[] {0, 1}),
				MatrixUtils.selectColumns(data2, new int[] {2, 3}),
				kNNs, expectedFromMILCA_2);
	}

	/**
	 * Test the computed multivariate MI against that calculated by Kraskov's own MILCA
	 * tool on the same data.
	 * 
	 * To run Kraskov's tool (http://www.klab.caltech.edu/~kraskov/MILCA/) for this 
	 * data, run:
	 * ./MIxnyn <dataFile> 2 2 3000 <kNearestNeighbours> 0
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
		double[] expectedFromMILCA_2 = {5.00322122, 4.29011291, 3.91312749, 
				3.69192886, 3.52807488, 3.39865354, 3.05327646, 2.79951639};
		
		System.out.println("Kraskov comparison 6 - multivariate dependent data 1");
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0, 1}),
				MatrixUtils.selectColumns(data, new int[] {2, 3}),
				kNNs, expectedFromMILCA_2);

	}

	/**
	 * Test the computed multivariate MI against that calculated by Kraskov's own MILCA
	 * tool on the same data.
	 * 
	 * To run Kraskov's tool (http://www.klab.caltech.edu/~kraskov/MILCA/) for this 
	 * data, run:
	 * ./MIxnyn <dataFile> 2 2 3000 <kNearestNeighbours> 0
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
		double[] expectedFromMILCA_2 = {0.33738970, 0.36251531, 0.34708687, 
				0.36200563, 0.35766125, 0.35007623, 0.35023664, 0.33728287};
		
		System.out.println("Kraskov comparison 7 - multivariate dependent data 1");
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0, 1}),
				MatrixUtils.selectColumns(data, new int[] {2, 3}),
				kNNs, expectedFromMILCA_2);

	}

	/**
	 * Test the computed multivariate MI against that calculated by Kraskov's own MILCA
	 * tool on the same data.
	 * 
	 * To run Kraskov's tool (http://www.klab.caltech.edu/~kraskov/MILCA/) for this 
	 * data, run:
	 * ./MIxnyn <dataFile> 5 5 10000 <kNearestNeighbours> 0
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
		double[] expectedFromMILCA_2 = {0.00815609, 0.00250864, 0.00035825,
				0.00172174, 0.00033354};
		
		System.out.println("Kraskov comparison 8 - multivariate uncorrelated Gaussian data 1");
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0, 1, 2, 3, 4}),
				MatrixUtils.selectColumns(data, new int[] {5, 6, 7, 8, 9}),
				kNNs, expectedFromMILCA_2);

	}

	/**
	 * Test the computed multivariate MI against that calculated by Kraskov's own MILCA
	 * tool on the same data.
	 * 
	 * To run Kraskov's tool (http://www.klab.caltech.edu/~kraskov/MILCA/) for this 
	 * data, run:
	 * ./MIxnyn <dataFile> 1 1 10000 <kNearestNeighbours> 0
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testMIforRandomGaussianVariablesFromLargeFile() throws Exception {
		
		// Test set 9:
		
		// We'll take the columns from this data set
		ArrayFileReader afr = new ArrayFileReader("demos/data/10ColsRandomGaussian-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {1, 2, 4, 10, 15};
		// Expected values from Kraskov's MILCA toolkit:
		double[] expectedFromMILCA_2 = {0.01542004, 0.01137151, 0.00210945,
				0.00159921, 0.00031277};
		
		System.out.println("Kraskov comparison 9 - uncorrelated Gaussian data 1 - large file");
		checkMIForGivenData(MatrixUtils.selectColumns(data, new int[] {0}),
				MatrixUtils.selectColumns(data, new int[] {1}),
				kNNs, expectedFromMILCA_2);

	}
}
