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
 * = Example 8 - Transfer entropy on continuous data by binning =
 * 
 * Simple transfer entropy (TE) calculation on continuous-valued data 
 * by binning the continuous data to discrete, then using a discrete TE calculator.
 * 
 * @author Joseph Lizier
 *
 */
public class Example8TeContinuousDataByBinning {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		
		// Prepare to generate some random normalised data.
		int numObservations = 1000;
		double covariance = 0.4;
		int numDiscreteLevels = 4;

		// Create destArray correlated to previous value of sourceArray:
		RandomGenerator rg = new RandomGenerator();
		double[] sourceArray = rg.generateNormalData(numObservations, 0, 1);
		double[] destArray = rg.generateNormalData(numObservations, 0, 1-covariance);
		for (int t = 1; t < numObservations; t++) {
			destArray[t] += covariance * sourceArray[t-1];
		}

		// Discretize or bin the data -- one could also call:
		//  MatrixUtils.discretiseMaxEntropy for a maximum entropy binning
		int[] binnedSource = MatrixUtils.discretise(sourceArray, numDiscreteLevels);
		int[] binnedDest = MatrixUtils.discretise(destArray, numDiscreteLevels);
		
		// Create a TE calculator and run it:
		TransferEntropyCalculatorDiscrete teCalc=
                new TransferEntropyCalculatorDiscrete(numDiscreteLevels, 1);
		teCalc.initialise();
		teCalc.addObservations(binnedSource, binnedDest);
		double result = teCalc.computeAverageLocalOfObservations();
		// Calculation will be heavily biased because of the binning,
		//  and the small number of samples
		System.out.printf("TE result %.4f bits; expected to be close to " +
				"%.4f bits for these correlated Gaussians\n",
				result, Math.log(1.0/(1-Math.pow(covariance,2)))/Math.log(2));
	}
}
