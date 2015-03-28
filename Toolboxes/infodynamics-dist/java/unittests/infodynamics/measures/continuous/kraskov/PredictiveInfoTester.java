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

import infodynamics.measures.continuous.PredictiveInfoAbstractTester;

public class PredictiveInfoTester extends PredictiveInfoAbstractTester {

	protected String NUM_THREADS_TO_USE_DEFAULT = MutualInfoCalculatorMultiVariateKraskov.USE_ALL_THREADS;
	protected String NUM_THREADS_TO_USE = NUM_THREADS_TO_USE_DEFAULT;
	
	/**
	 * Utility function to create a calculator for the given algorithm number
	 * 
	 * @param algNumber
	 * @return
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 */
	public PredictiveInfoCalculatorKraskov getNewCalc(int algNumber)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		return new PredictiveInfoCalculatorKraskov(algNumber);
	}
	
	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 * 
	 */
	public void checkLocalsAverageCorrectly(int algNumber, String numThreads) throws Exception {
		
		PredictiveInfoCalculatorKraskov piCalc = getNewCalc(algNumber);
		
		String kraskov_K = "4";
		
		piCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);
		piCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
				numThreads);

		super.testLocalsAverageCorrectly(piCalc, 2, 10000);
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
		
		PredictiveInfoCalculatorKraskov piCalc = getNewCalc(algNumber);
		
		String kraskov_K = "4";
		
		piCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_K,
				kraskov_K);
		piCalc.setProperty(
				MutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
				NUM_THREADS_TO_USE);

		super.testComputeSignificanceDoesntAlterAverage(piCalc, 2, 100);
	}
	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		checkComputeSignificanceDoesntAlterAverage(1);
		checkComputeSignificanceDoesntAlterAverage(2);
	}
}
