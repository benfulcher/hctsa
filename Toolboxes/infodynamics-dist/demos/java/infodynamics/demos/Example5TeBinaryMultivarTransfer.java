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

package infodynamics.demos;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;
import infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete;

/**
 * 
 * = Example 5 - Multivariate transfer entropy on binary data =
 * 
 * Multivariate transfer entropy (TE) calculation on binary data using the discrete TE calculator:
 * 
 * @author Joseph Lizier
 *
 */
public class Example5TeBinaryMultivarTransfer {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		
		// Generate some random binary data.
		int timeSeriesLength = 100;
		RandomGenerator rg = new RandomGenerator();
		int[][] sourceArray = rg.generateRandomInts(timeSeriesLength, 2, 2);
		int[][] sourceArray2 = rg.generateRandomInts(timeSeriesLength, 2, 2);
		// Destination variable takes a copy of the first bit of the 
		//  previous source value in bit 0,
		//  and an XOR of the two previous bits of the source in bit 1:
		int[][] destArray = new int[timeSeriesLength][2];
		for (int r = 1; r < timeSeriesLength; r++) {
			// This is a bitwise XOR, but is fine for our purposes
			//  with binary data:
			destArray[r][0] = sourceArray[r - 1][0];
			destArray[r][1]	= sourceArray[r - 1][0] ^ sourceArray[r - 1][1];
		}

		// Create a TE calculator and run it.
		// Need to represent 4-state variables for the joint destination variable
		TransferEntropyCalculatorDiscrete teCalc=
				new TransferEntropyCalculatorDiscrete(4, 1);
		teCalc.initialise();

		// We need to construct the joint values of the dest and source before we pass them in:
		teCalc.addObservations(MatrixUtils.computeCombinedValues(sourceArray, 2),
				MatrixUtils.computeCombinedValues(destArray, 2));

		double result = teCalc.computeAverageLocalOfObservations();
		System.out.printf("For source which the 2 bits are determined from, " +
				"result should be close to 2 bits : %.3f\n", result);
		
		// Check random source:
		teCalc.initialise();
		teCalc.addObservations(MatrixUtils.computeCombinedValues(sourceArray2, 2),
				MatrixUtils.computeCombinedValues(destArray, 2));
		double result2 = teCalc.computeAverageLocalOfObservations();
		System.out.printf("For random source, result should be close to 0 bits " +
				"in theory: %.3f\n", result2);
		
		System.out.printf("The result for random source is inflated towards 0.3 " +
				" due to finite observation length (%d). One can verify that the " +
				"answer is consistent with that from a random source by checking: " +
				"teCalc.computeSignificance(1000); ans.pValue\n",
				teCalc.getNumObservations());
	}
}
