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
import infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov;

/**
 * 
 * = Example 4 - Transfer entropy on continuous data using Kraskov estimators =
 * 
 * Simple transfer entropy (TE) calculation on continuous-valued data using the Kraskov-estimator TE calculator.
 * 
 * @author Joseph Lizier
 *
 */
public class Example4TeContinuousDataKraskov {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		
		// Generate some random normalised data.
		int numObservations = 1000;
		double covariance = 0.4;

		// Create destArray correlated to previous value of sourceArray:
		RandomGenerator rg = new RandomGenerator();
		double[] sourceArray = rg.generateNormalData(numObservations, 0, 1);
		double[] destArray = rg.generateNormalData(numObservations, 0, 1-covariance);
		for (int t = 1; t < numObservations; t++) {
			destArray[t] += covariance * sourceArray[t-1];
		}
		// And an uncorrelated second source
		double[] sourceArray2 = rg.generateNormalData(numObservations, 0, 1);

		// Create a TE calculator and run it:
		TransferEntropyCalculatorKraskov teCalc =
				new TransferEntropyCalculatorKraskov();
		teCalc.setProperty("k", "4"); // Use Kraskov parameter K=4 for 4 nearest neighbours
		teCalc.initialise(1); // Use history length 1 (Schreiber k=1)
		
		// Perform calculation with correlated source:
		teCalc.setObservations(sourceArray, destArray);
		double result = teCalc.computeAverageLocalOfObservations();
		// Note that the calculation is a random variable (because the generated
		//  data is a set of random variables) - the result will be of the order
		//  of what we expect, but not exactly equal to it; in fact, there will
		//  be a large variance around it.
		System.out.printf("TE result %.4f nats; expected to be close to " +
				"%.4f nats for these correlated Gaussians\n",
				result, Math.log(1.0/(1-Math.pow(covariance,2))));

		//  Perform calculation with uncorrelated source:
		teCalc.initialise(); // Initialise leaving the parameters the same
		teCalc.setObservations(sourceArray2, destArray);
		// For random source, it should give something close to 0 bits
		double result2 = teCalc.computeAverageLocalOfObservations();
		System.out.printf("TE result %.4f nats; expected to be close to " +
				"0 nats for these uncorrelated Gaussians\n", result2);
		
		// We can also compute the local TE values for the time-series samples here:
		//  (See more about utility of local TE in the CA demos)
		double[] localTE = teCalc.computeLocalOfPreviousObservations();
		System.out.printf("Notice that the mean of locals, %.4f nats," +
				" equals the previous result\n",
				MatrixUtils.sum(localTE)/(double)(numObservations-1));
	}
}
