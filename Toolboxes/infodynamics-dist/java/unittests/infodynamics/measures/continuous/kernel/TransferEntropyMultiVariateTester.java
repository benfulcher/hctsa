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

package infodynamics.measures.continuous.kernel;

public class TransferEntropyMultiVariateTester
	extends infodynamics.measures.continuous.TransferEntropyMultiVariateAbstractTester {

	/**
	 * Confirm that the local values average correctly back to the average value
	 * 
	 */
	public void testLocalsAverageCorrectly() throws Exception {
		
		TransferEntropyCalculatorMultiVariateKernel teCalc =
				new TransferEntropyCalculatorMultiVariateKernel();
		
		String kernelWidth = "1";
		
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKernel.NORMALISE_PROP_NAME,
				"true");
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKernel.EPSILON_PROP_NAME,
				kernelWidth);

		super.testLocalsAverageCorrectly(teCalc, 2, 100, 1);
	}
	
	/**
	 * Confirm that significance testing doesn't alter the average that
	 * would be returned.
	 * 
	 * @throws Exception
	 */
	public void testComputeSignificanceDoesntAlterAverage() throws Exception {
		
		TransferEntropyCalculatorMultiVariateKernel teCalc =
				new TransferEntropyCalculatorMultiVariateKernel();
		
		String kernelWidth = "1";
		
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKernel.NORMALISE_PROP_NAME,
				"true");
		teCalc.setProperty(
				TransferEntropyCalculatorMultiVariateKernel.EPSILON_PROP_NAME,
				kernelWidth);

		super.testComputeSignificanceDoesntAlterAverage(teCalc, 2, 100, 1);
	}

}
