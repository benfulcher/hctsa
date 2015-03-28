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

package infodynamics.measures.continuous.gaussian;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;

public class EntropyCalculatorGaussianTest extends TestCase {

	public void testVarianceSetting() {
		EntropyCalculatorGaussian entcalc = new EntropyCalculatorGaussian();
		entcalc.initialise();
		double variance = 3.45567;
		entcalc.setVariance(variance);
		assertEquals(variance, entcalc.variance);
	}

	public void testVarianceCalculation() {
		EntropyCalculatorGaussian entcalc = new EntropyCalculatorGaussian();
		entcalc.initialise();
		RandomGenerator randomGenerator = new RandomGenerator();
		double[] gaussianObservations = randomGenerator.generateNormalData(100, 5, 3);
		double variance = MatrixUtils.stdDev(gaussianObservations);
		variance *= variance;
		entcalc.setObservations(gaussianObservations);
		assertEquals(variance, entcalc.variance);
	}
	
	public void testEntropyCalculation() {
		EntropyCalculatorGaussian entcalc = new EntropyCalculatorGaussian();
		entcalc.initialise();
		double variance = 3.45567;
		entcalc.setVariance(variance);
		double expectedEntropy = 0.5 * Math.log(2.0*Math.PI*Math.E*variance);
		assertEquals(expectedEntropy, entcalc.computeAverageLocalOfObservations());
	}
}
