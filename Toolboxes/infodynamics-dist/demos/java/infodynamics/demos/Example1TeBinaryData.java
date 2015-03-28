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

import infodynamics.utils.RandomGenerator;
import infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete;

/**
 * 
 * = Example 1 - Transfer entropy on binary data =
 * 
 * Simple transfer entropy (TE) calculation on binary data using the discrete TE calculator.
 * 
 * @author Joseph Lizier
 *
 */
public class Example1TeBinaryData {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		int arrayLengths = 100;
		RandomGenerator rg = new RandomGenerator();
		
		// Generate some random binary data:
		int[] sourceArray = rg.generateRandomInts(arrayLengths, 2);
		int[] destArray = new int[arrayLengths];
		destArray[0] = 0;
		System.arraycopy(sourceArray, 0, destArray, 1, arrayLengths - 1);
		int[] sourceArray2 = rg.generateRandomInts(arrayLengths, 2);
		
		// Create a TE calculator and run it:
		TransferEntropyCalculatorDiscrete teCalc=
				new TransferEntropyCalculatorDiscrete(2, 1);
		teCalc.initialise();
		teCalc.addObservations(sourceArray, destArray);
		double result = teCalc.computeAverageLocalOfObservations();
		System.out.printf("For copied source, result should be close to 1 bit : %.3f bits\n", result);
		teCalc.initialise();
		teCalc.addObservations(sourceArray2, destArray);
		double result2 = teCalc.computeAverageLocalOfObservations();
		System.out.printf("For random source, result should be close to 0 bits: %.3f bits\n", result2);
	}

}
