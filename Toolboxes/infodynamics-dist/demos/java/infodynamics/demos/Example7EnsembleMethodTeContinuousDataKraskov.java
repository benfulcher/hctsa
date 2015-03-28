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
import infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov;

/**
 * 
 * = Example 7 - Ensemble method with transfer entropy on continuous data using Kraskov estimators =
 * 
 * Calculation of transfer entropy (TE) by supplying an ensemble of
 *  samples from multiple time series.
 * We use continuous-valued data using the Kraskov-estimator TE calculator
 *  here.
 * 
 * @author Joseph Lizier
 *
 */
public class Example7EnsembleMethodTeContinuousDataKraskov {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		
		// Prepare to generate some random normalised data.
		int numObservations = 1000;
		double covariance = 0.4;
		RandomGenerator rg = new RandomGenerator();

		// Create a TE calculator and run it:
		TransferEntropyCalculatorKraskov teCalc =
				new TransferEntropyCalculatorKraskov();
		teCalc.setProperty("k", "4"); // Use Kraskov parameter K=4 for 4 nearest neighbours
		teCalc.initialise(1); // Use history length 1 (Schreiber k=1)
		teCalc.startAddObservations();
		
		for (int trial = 0; trial < 10; trial++) {
		
			// Create a new trial, with destArray correlated to
			//  previous value of sourceArray:
			double[] sourceArray = rg.generateNormalData(numObservations, 0, 1);
			double[] destArray = rg.generateNormalData(numObservations, 0, 1-covariance);
			for (int t = 1; t < numObservations; t++) {
				destArray[t] += covariance * sourceArray[t-1];
			}
			
			// Add observations for this trial:
			System.out.printf("Adding samples from trial %d ...\n", trial);
			teCalc.addObservations(sourceArray, destArray);
		}
		
		// We've finished adding trials:
		System.out.println("Finished adding trials");
		teCalc.finaliseAddObservations();
		
		// Compute the result:
		System.out.println("Computing TE ...");
		double result = teCalc.computeAverageLocalOfObservations();
		// Note that the calculation is a random variable (because the generated
		//  data is a set of random variables) - the result will be of the order
		//  of what we expect, but not exactly equal to it; in fact, there will
		//  be some variance around it (smaller than example 4 since we have more samples).
		System.out.printf("TE result %.4f nats; expected to be close to " +
				"%.4f nats for these correlated Gaussians\n",
				result, Math.log(1.0/(1-Math.pow(covariance,2))));
	}
}
