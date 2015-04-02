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
 * = Example 2 - Transfer entropy on multidimensional binary data =
 * 
 * Simple transfer entropy (TE) calculation on multidimensional binary data using the discrete TE calculator.
 *
 * This example shows how to handle multidimensional arrays where
 *  we pool the observations over all variables with the discrete calculator.
 * 
 * @author Joseph Lizier
 *
 */
public class Example2TeMultidimBinaryData {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		int timeSteps = 2;
		int variables = 100;
		RandomGenerator rg = new RandomGenerator();
		
		// Create many columns in a multidimensional array (2 rows by 100 columns),
		//  where the next time step (row 2) copies the value of the column on the left
		//  from the previous time step (row 1):
		int[][] twoDTimeSeries = new int[timeSteps][];
		twoDTimeSeries[0] = rg.generateRandomInts(variables, 2);
		twoDTimeSeries[1] = new int[variables];
		twoDTimeSeries[1][0] = twoDTimeSeries[0][variables - 1];
		System.arraycopy(twoDTimeSeries[0], 0, twoDTimeSeries[1], 1, variables - 1);

		// Create a TE calculator and run it:
		TransferEntropyCalculatorDiscrete teCalc=
				new TransferEntropyCalculatorDiscrete(2, 1);
		teCalc.initialise();
		// Add observations of transfer across one cell to the right (j=1)
		//  per time step:
		teCalc.addObservations(twoDTimeSeries, 1);

		double result2D = teCalc.computeAverageLocalOfObservations();
		System.out.printf("The result should be close to 1 bit here, " +
				"since we are executing copy operations of what is effectively " +
				"a random bit to each cell here: %.3f bits\n", result2D);
	}
}
