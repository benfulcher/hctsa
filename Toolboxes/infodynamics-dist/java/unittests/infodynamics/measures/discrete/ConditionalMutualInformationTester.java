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

package infodynamics.measures.discrete;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.RandomGenerator;

import junit.framework.TestCase;

public class ConditionalMutualInformationTester extends TestCase {

	protected RandomGenerator rand = new RandomGenerator();
	
	protected int numObservations = 1000;
	
	public void testComputeSignificanceInt() {
		// Just making sure that no exception is thrown
		ConditionalMutualInformationCalculatorDiscrete condMiCalc = new ConditionalMutualInformationCalculatorDiscrete(2, 2, 2);
		int[] x1 = rand.generateRandomInts(numObservations, 2);
		int[] x2 = rand.generateRandomInts(numObservations, 2);
		int[] cond = rand.generateRandomInts(numObservations, 2);
		condMiCalc.initialise();
		condMiCalc.addObservations(x1, x2, cond);
		condMiCalc.computeAverageLocalOfObservations();
		condMiCalc.computeSignificance(1000);
	}

	public void testComputeLocalsGivesCorrectAverage() {
		ConditionalMutualInformationCalculatorDiscrete condMiCalc = new ConditionalMutualInformationCalculatorDiscrete(2, 2, 2);
		int[] x1 = rand.generateRandomInts(numObservations, 2);
		int[] x2 = rand.generateRandomInts(numObservations, 2);
		int[] cond = rand.generateRandomInts(numObservations, 2);
		condMiCalc.initialise();
		condMiCalc.addObservations(x1, x2, cond);
		double average = condMiCalc.computeAverageLocalOfObservations();
		double avLastAverage = condMiCalc.getLastAverage();
		double[] locals = condMiCalc.computeLocal(x1, x2, cond);
		double localsLastAverage = condMiCalc.getLastAverage();
		assertEquals(average, MatrixUtils.mean(locals), 0.000001);
		assertEquals(avLastAverage, localsLastAverage, 0.000001);
		assertEquals(average, avLastAverage, 0.000001);
	}

	public void testSetDebug() {
		ConditionalMutualInformationCalculatorDiscrete condMiCalc = new ConditionalMutualInformationCalculatorDiscrete(2, 2, 2);
		assertFalse(condMiCalc.debug);
		condMiCalc.setDebug(true);
		assertTrue(condMiCalc.debug);
		condMiCalc.setDebug(false);
		assertFalse(condMiCalc.debug);
	}

}
