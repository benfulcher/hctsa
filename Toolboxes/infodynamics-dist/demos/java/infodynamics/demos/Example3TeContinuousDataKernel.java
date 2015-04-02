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

import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;
import infodynamics.measures.continuous.kernel.TransferEntropyCalculatorKernel;

/**
 * 
 * = Example 3 - Transfer entropy on continuous data using kernel estimators =
 * 
 * Simple transfer entropy (TE) calculation on continuous-valued data using the (box) kernel-estimator TE calculator.
 * 
 * @author Joseph Lizier
 *
 */
public class Example3TeContinuousDataKernel {

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
		TransferEntropyCalculatorKernel teCalc =
				new TransferEntropyCalculatorKernel();
		teCalc.setProperty("NORMALISE", "true"); // Normalise the individual variables (default)
		teCalc.initialise(1, 0.5); // Use history length 1 (Schreiber k=1), kernel width of 0.5 normalised units
		teCalc.setObservations(sourceArray, destArray);
		// For copied source, should give something close to expected value for correlated Gaussians:
		double result = teCalc.computeAverageLocalOfObservations();
		System.out.printf("TE result %.4f bits; expected to be close to " +
				"%.4f bits for these correlated Gaussians but biased upwards\n",
				result, Math.log(1.0/(1-Math.pow(covariance,2)))/Math.log(2));

		teCalc.initialise(); // Initialise leaving the parameters the same
		teCalc.setObservations(sourceArray2, destArray);
		// For random source, it should give something close to 0 bits
		double result2 = teCalc.computeAverageLocalOfObservations();
		System.out.printf("TE result %.4f bits; expected to be close to " +
				"0 bits for uncorrelated Gaussians but will be biased upwards\n",
				result2);
		
		// We can get insight into the bias by examining the null distribution:
		EmpiricalMeasurementDistribution nullDist = teCalc.computeSignificance(100);
		System.out.printf("Null distribution for unrelated source and destination " +
				"(i.e. the bias) has mean %.4f and standard deviation %.4f\n",
				nullDist.getMeanOfDistribution(), nullDist.getStdOfDistribution());
	}
}
