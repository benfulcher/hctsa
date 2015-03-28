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

public class TransferEntropyMultiVariateTester
	extends infodynamics.measures.continuous.TransferEntropyMultiVariateAbstractTester {

	protected String NUM_THREADS_TO_USE_DEFAULT = ConditionalMutualInfoCalculatorMultiVariateKraskov.USE_ALL_THREADS;
	protected String NUM_THREADS_TO_USE = NUM_THREADS_TO_USE_DEFAULT;
	
	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 */
	public void testLocalsAverageCorrectly() throws Exception {
		
		TransferEntropyCalculatorMultiVariateKraskov teCalc =
				new TransferEntropyCalculatorMultiVariateKraskov();
		
		String kraskov_K = "4";
		
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKraskov.PROP_KRASKOV_ALG_NUM,
				"2");
		teCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);

		super.testLocalsAverageCorrectly(teCalc, 2, 100, 1);
	}
	
	/**
	 * Confirm that significance testing doesn't alter the average that
	 * would be returned.
	 * 
	 * @throws Exception
	 */
	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		
		TransferEntropyCalculatorMultiVariateKraskov teCalc =
				new TransferEntropyCalculatorMultiVariateKraskov();
		
		String kraskov_K = "4";
		
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKraskov.PROP_KRASKOV_ALG_NUM,
				"2");
		teCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);

		super.testComputeSignificanceDoesntAlterAverage(teCalc, 2, 100, 1);
	}

	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 */
	public void testUnivariateSignatureMatchesMultivariate() throws Exception {
		
		TransferEntropyCalculatorMultiVariateKraskov teCalc =
				new TransferEntropyCalculatorMultiVariateKraskov();
		
		String kraskov_K = "4";
		
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKraskov.PROP_KRASKOV_ALG_NUM,
				"1");
		teCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);

		super.testUnivariateMatchesMultivariateRoute(teCalc, 100, 1);
	}

	/**
	 * Utility function to create a calculator for the given algorithm number
	 * 
	 * @param algNumber
	 * @return
	 */
	public TransferEntropyCalculatorKraskov getNewCalc(int algNumber) throws Exception {
		TransferEntropyCalculatorKraskov teCalc =
				new TransferEntropyCalculatorKraskov();
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_KRASKOV_ALG_NUM,
					Integer.toString(algNumber));
		return teCalc;
	}
	
	/**
	 * Utility function to run Kraskov algorithm 1
	 *  as transfer entropy for data with known results
	 *  from TRENTOOL. (with default parameter settings k=1, l=1)
	 * 
	 * @param var1 source time-series data set
	 * @param var2 dest time-series data set
	 * @param kNNs array of Kraskov k nearest neighbours parameter to check
	 * @param expectedResults array of expected results for each k
	 */
	protected void checkTEForGivenData(double[] var1, double[] var2,
			int[] kNNs, double[] expectedResults) throws Exception {
		checkTEForGivenData(var1, var2, 1, 1, kNNs, expectedResults);
	}
	
	/**
	 * Utility function to run Kraskov algorithm 1
	 *  as transfer entropy for data with known results
	 *  from TRENTOOL.
	 * 
	 * @param var1 source time-series data set
	 * @param var2 dest time-series data set
	 * @param historyK history length k of destination
	 * @param historyL history length l of source
	 * @param kNNs array of Kraskov k nearest neighbours parameter to check
	 * @param expectedResults array of expected results for each k
	 */
	protected void checkTEForGivenData(double[] var1, double[] var2,
			int historyK, int historyL, int[] kNNs, double[] expectedResults) throws Exception {
				
		TransferEntropyCalculatorKraskov teCalc = getNewCalc(1);
		
		// Normalise the data ourselves rather than letting the calculator do it -
		//  this ensures the extra values in the time series (e.g. last value in source)
		//  are taken into account, in line with TRENTOOL
		var1 = MatrixUtils.normaliseIntoNewArray(var1);
		var2 = MatrixUtils.normaliseIntoNewArray(var2);
		
		for (int kIndex = 0; kIndex < kNNs.length; kIndex++) {
			int k = kNNs[kIndex];
			teCalc.setProperty(
					ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K,
					Integer.toString(k));
			// We already normalised above, and this will do a different
			//  normalisation without taking the extra values in to account if we did it
			teCalc.setProperty(
					ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORMALISE,
					Boolean.toString(false));
			// No longer need to set this property as it's set by default:
			//teCalc.setProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE,
			//		EuclideanUtils.NORM_MAX_NORM_STRING);
			teCalc.setProperty(
					MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
					NUM_THREADS_TO_USE);
			teCalc.initialise(historyK, 1, historyL, 1, 1);
			// And set the observations
			teCalc.setObservations(var1, var2);
			double te = teCalc.computeAverageLocalOfObservations();
			//teCalc.setDebug(false);
			
			System.out.printf("k=%d: Average TE %.8f (expected %.8f)\n",
					k, te, expectedResults[kIndex]);
			// 6 decimal places is Matlab accuracy
			assertEquals(expectedResults[kIndex], te, 0.000001);			
		}
	}

	/**
	 * Test the computed univariate TE
	 * against that calculated by Wibral et al.'s TRENTOOL
	 * on the same data.
	 * 
	 * To run TRENTOOL (http://www.trentool.de/) for this 
	 * data, run its TEvalues.m matlab script on the multivariate source
	 * and dest data sets as:
	 * TEvalues(source, dest, 1, 1, 1, kraskovK, 0)
	 * with these values ensuring source-dest lag 1, history k=1,
	 * embedding lag 1, no dynamic correlation exclusion 
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testUnivariateTEforCoupledVariablesFromFile() throws Exception {
		
		// Test set 1:
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/2coupledRandomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {4};
		// Expected values from TRENTOOL:
		double[] expectedFromTRENTOOL = {0.3058006};
		
		System.out.println("Kraskov TE comparison 1 - univariate coupled data 1");
		checkTEForGivenData(MatrixUtils.selectColumn(data, 0),
				MatrixUtils.selectColumn(data, 1),
				kNNs, expectedFromTRENTOOL);
		
		// And now in the reverse direction:
		expectedFromTRENTOOL = new double[] {-0.0029744};
		
		System.out.println("  reverse direction:");
		checkTEForGivenData(MatrixUtils.selectColumn(data, 1),
				MatrixUtils.selectColumn(data, 0),
				kNNs, expectedFromTRENTOOL);

	}

	/**
	 * Test the computed univariate TE
	 * against that calculated by Wibral et al.'s TRENTOOL
	 * on the same data.
	 * 
	 * To run TRENTOOL (http://www.trentool.de/) for this 
	 * data, run its TEvalues.m matlab script on the multivariate source
	 * and dest data sets as:
	 * TEvalues(source, dest, 1, 1, 1, kraskovK, 0)
	 * with these values ensuring source-dest lag 1, history k=1,
	 * embedding lag 1, no dynamic correlation exclusion 
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testUnivariateTEforRandomDataFromFile() throws Exception {
		
		// Test set 2:
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/4randomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {4};
		// Expected values from TRENTOOL:
		double[] expectedFromTRENTOOL = {-0.0096556};
		
		System.out.println("Kraskov TE comparison 2 - univariate random data 1 (col 0->1)");
		checkTEForGivenData(MatrixUtils.selectColumn(data, 0),
				MatrixUtils.selectColumn(data, 1),
				kNNs, expectedFromTRENTOOL);
		
		// And now for other columns
		expectedFromTRENTOOL = new double[] {0.0175389};
		
		System.out.println("  (col 1->2):");
		checkTEForGivenData(MatrixUtils.selectColumn(data, 1),
				MatrixUtils.selectColumn(data, 2),
				kNNs, expectedFromTRENTOOL);

		// And now for other columns
		expectedFromTRENTOOL = new double[] {0.0026367};
		
		System.out.println("  (col 1->0):");
		checkTEForGivenData(MatrixUtils.selectColumn(data, 1),
				MatrixUtils.selectColumn(data, 0),
				kNNs, expectedFromTRENTOOL);

		// And now for other columns
		expectedFromTRENTOOL = new double[] {-0.00012474};
		
		System.out.println("  (col 0->2):");
		checkTEForGivenData(MatrixUtils.selectColumn(data, 0),
				MatrixUtils.selectColumn(data, 2),
				kNNs, expectedFromTRENTOOL);

		// And now for other columns
		expectedFromTRENTOOL = new double[] {-5.4437e-03};
		
		System.out.println("  (col 2->0):");
		checkTEForGivenData(MatrixUtils.selectColumn(data, 2),
				MatrixUtils.selectColumn(data, 0),
				kNNs, expectedFromTRENTOOL);
	}

	/**
	 * Test the computed multivariate conditional MI
	 * against that calculated by Wibral et al.'s TRENTOOL
	 * on the same data.
	 * 
	 * It's multivariate because we use embedding dimension 2 on both source
	 *  and destination.
	 * 
	 * To run TRENTOOL (http://www.trentool.de/) for this 
	 * data, run its TEvalues.m matlab script on the multivariate source
	 * and dest data sets as:
	 * TEvalues(source, dest, 2, 2, 1, kraskovK, 0)
	 * with these values ensuring source-dest lag 1, history k=2,
	 * history embedding dimension l=2 on source as well.
	 * embedding lag 1, no dynamic correlation exclusion 
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testMultivariateCondMIforCoupledDataFromFile() throws Exception {
		
		// Test set 3:
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/4ColsPairedOneStepNoisyDependence-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {4};
		// Expected values from TRENTOOL:
		double[] expectedFromTRENTOOL = {0.1400645};
		
		System.out.println("Kraskov TE - multivariate coupled data 1, k=2,l=2");
		System.out.println("  (0->2)");
		checkTEForGivenData(MatrixUtils.selectColumn(data, 0),
				MatrixUtils.selectColumn(data, 2),
				2, 2,
				kNNs, expectedFromTRENTOOL);
		
		// And now for reverse direction:
		expectedFromTRENTOOL = new double[] {-0.0181459};
		
		System.out.println("  (2->0):");
		checkTEForGivenData(MatrixUtils.selectColumn(data, 2),
				MatrixUtils.selectColumn(data, 0),
				2, 2,
				kNNs, expectedFromTRENTOOL);		

		// And now for other columns:
		expectedFromTRENTOOL = new double[] {0.1639186};
		System.out.println("  (1->3):");
		checkTEForGivenData(MatrixUtils.selectColumn(data, 1),
				MatrixUtils.selectColumn(data, 3),
				2, 2,
				kNNs, expectedFromTRENTOOL);		
		// And in reverse:
		expectedFromTRENTOOL = new double[] {0.0036976};
		System.out.println("  (3->1):");
		checkTEForGivenData(MatrixUtils.selectColumn(data, 3),
				MatrixUtils.selectColumn(data, 1),
				2, 2,
				kNNs, expectedFromTRENTOOL);
		
		// -------------
		// And finally, confirm that we get different results for k=1,l=1,
		//  which match TRENTOOL
		expectedFromTRENTOOL = new double[] {0.0072169};
		System.out.println("  (0->1) but with k=1,l=1:");
		checkTEForGivenData(MatrixUtils.selectColumn(data, 0),
				MatrixUtils.selectColumn(data, 1),
				1, 1,
				kNNs, expectedFromTRENTOOL);		
		// And in reverse
		expectedFromTRENTOOL = new double[] {0.0011738};
		System.out.println("  (1->2) but with k=1,l=1:");
		checkTEForGivenData(MatrixUtils.selectColumn(data, 1),
				MatrixUtils.selectColumn(data, 2),
				1, 1,
				kNNs, expectedFromTRENTOOL);		
	}

	/**
	 * Test the computed univariate TE
	 * using various numbers of threads.
	 * 
	 * Test the computed univariate TE
	 * against that calculated by Wibral et al.'s TRENTOOL
	 * on the same data.
	 * 
	 * To run TRENTOOL (http://www.trentool.de/) for this 
	 * data, run its TEvalues.m matlab script on the multivariate source
	 * and dest data sets as:
	 * TEvalues(source, dest, 1, 1, 1, kraskovK, 0)
	 * with these values ensuring source-dest lag 1, history k=1,
	 * embedding lag 1, no dynamic correlation exclusion 
	 * 
	 * @throws Exception if file not found 
	 * 
	 */
	public void testUnivariateTEVariousNumberThreads() throws Exception {
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/4randomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();
		
		// Use various Kraskov k nearest neighbours parameter
		int[] kNNs = {4};
		// Expected values from TRENTOOL:
		double[] expectedFromTRENTOOL = {-0.0096556};
		
		System.out.println("Kraskov TE - multivariate coupled data 1, k=2,l=2 (0->2)");
		System.out.println(" with various numbers of threads:");
		System.out.println(" -- 1 thread:");
		NUM_THREADS_TO_USE = "1";
		checkTEForGivenData(MatrixUtils.selectColumn(data, 0),
				MatrixUtils.selectColumn(data, 1),
				kNNs, expectedFromTRENTOOL);
		System.out.println(" -- 2 threads:");
		NUM_THREADS_TO_USE = "2";
		checkTEForGivenData(MatrixUtils.selectColumn(data, 0),
				MatrixUtils.selectColumn(data, 1),
				kNNs, expectedFromTRENTOOL);
		System.out.println(" -- 3 threads:");
		NUM_THREADS_TO_USE = "3";
		checkTEForGivenData(MatrixUtils.selectColumn(data, 0),
				MatrixUtils.selectColumn(data, 1),
				kNNs, expectedFromTRENTOOL);
		System.out.println(" -- all threads:");
		NUM_THREADS_TO_USE = ConditionalMutualInfoCalculatorMultiVariateKraskov.USE_ALL_THREADS;
		checkTEForGivenData(MatrixUtils.selectColumn(data, 0),
				MatrixUtils.selectColumn(data, 1),
				kNNs, expectedFromTRENTOOL);
		
		// And finally test that multithreading is still ok if we have
		//  an imbalanced number of data between each thread.
		// Expected value is only generated from our own code; we're not so 
		//  interested in checking for this precise value, as we are in
		//  checking that the value is stable when we change the number of
		//  threads and have an uneven amount of data in each thread.
		double[] expectedValue = new double[] {0.026517704};
		NUM_THREADS_TO_USE = "2";
		checkTEForGivenData(
				MatrixUtils.select(
						MatrixUtils.selectColumn(data, 0),
						0, 501),
				MatrixUtils.select(
						MatrixUtils.selectColumn(data, 1),
						0, 501),
				kNNs, expectedValue);
		NUM_THREADS_TO_USE = "3";
		checkTEForGivenData(
				MatrixUtils.select(
						MatrixUtils.selectColumn(data, 0),
						0, 501),
				MatrixUtils.select(
						MatrixUtils.selectColumn(data, 1),
						0, 501),
				kNNs, expectedValue);
		NUM_THREADS_TO_USE = NUM_THREADS_TO_USE_DEFAULT;
	}
}
