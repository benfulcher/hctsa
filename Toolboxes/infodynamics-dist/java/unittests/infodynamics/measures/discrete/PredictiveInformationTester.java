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

import infodynamics.utils.RandomGenerator;

import junit.framework.TestCase;
import java.util.Random;

public class PredictiveInformationTester extends TestCase {

	public void testFullyDependent() {
		PredictiveInformationCalculatorDiscrete piCalc = new PredictiveInformationCalculatorDiscrete(2, 1);
		
		// Next row is the inverse of the one above
		int[] x = new int[101];
		for (int t = 1; t < 101; t++) {
			x[t] = (x[t-1] == 1) ? 0 : 1;
		}
		piCalc.initialise();
		piCalc.addObservations(x);
		double piInverses = piCalc.computeAverageLocalOfObservations();
		assertEquals(1.0, piInverses, 0.000000001);
	}

	public void testNoActivity() {
		PredictiveInformationCalculatorDiscrete piCalc = new PredictiveInformationCalculatorDiscrete(2, 1);
		
		int[] x = new int[101];
		piCalc.initialise();
		piCalc.addObservations(x);
		double piNoActivity = piCalc.computeAverageLocalOfObservations();
		assertEquals(0.0, piNoActivity, 0.000000001);		
	}

	public void testReinitialisation() {
		PredictiveInformationCalculatorDiscrete piCalc = new PredictiveInformationCalculatorDiscrete(2, 1);
		
		int[] timeSteps = new int[] {11, 101, 1001, 10001};
		
		for (int tsIndex = 0; tsIndex < timeSteps.length; tsIndex++) {
			// Next row is the inverse of the one above
			int[] x = new int[timeSteps[tsIndex]];
			for (int t = 1; t < timeSteps[tsIndex]; t++) {
				x[t] = (x[t-1] == 1) ? 0 : 1;
			}
			piCalc.initialise();
			piCalc.addObservations(x);
			double piInverses = piCalc.computeAverageLocalOfObservations();
			assertEquals(1.0, piInverses, 0.000000001);
			// Now there is no activity on x
			x = new int[timeSteps[tsIndex]];
			piCalc.initialise();
			piCalc.addObservations(x);
			double piNoActivity = piCalc.computeAverageLocalOfObservations();
			assertEquals(0.0, piNoActivity, 0.000000001);
		}
	}

	public void testIndependent() {
		PredictiveInformationCalculatorDiscrete piCalc = new PredictiveInformationCalculatorDiscrete(2, 1);
		
		// Next value is the independent of the previous
		int[] x = new int[] {0, 0, 1, 1, 0};
		piCalc.initialise();
		piCalc.addObservations(x);
		double piIndpt = piCalc.computeAverageLocalOfObservations();
		assertEquals(0, piIndpt, 0.000000001);
	}

	public void testConvergenceWithActiveInfoStorage() {
		RandomGenerator rg = new RandomGenerator();
		Random random = new Random();
		
		int[][] x = new int[100][100];
		// Initialise first row
		x[0] = rg.generateRandomInts(100, 2);
		for (int t = 1; t < 100; t++) {
			for (int c = 0; c < 100; c++) {
				// Copy the previous bit with some chance, else
				//  assign at random. This ensures some non-zero 
				//  active info storage
				x[t][c] = (Math.random() < 0.5) ? x[t-1][c] : random.nextInt(2);
			}
		}
		// Compute the predictive information and check that it
		//  matches the active info storage when both are calculated
		//  with history length 1.
		PredictiveInformationCalculatorDiscrete piCalc = new PredictiveInformationCalculatorDiscrete(2, 1);
		piCalc.initialise();
		piCalc.addObservations(x);
		double pi = piCalc.computeAverageLocalOfObservations();
		ActiveInformationCalculatorDiscrete aiCalc = new ActiveInformationCalculatorDiscrete(2, 1);
		aiCalc.initialise();
		aiCalc.addObservations(x);
		double ai = aiCalc.computeAverageLocalOfObservations();
		assertEquals(ai, pi, 0.000000001);
		System.out.printf("PI: %.5f == AI: %.5f\n", pi, ai);
	}

	public void testDetectionOfLongerTermTrends() {
		int[][] x = new int[100][4];
		// Initialise first two rows
		x[0] = new int[] {0, 0, 1, 1};
		x[1] = new int[] {0, 1, 0, 1};
		for (int t = 2; t < 100; t++) {
			for (int c = 0; c < 4; c++) {
				// Copy the bit two steps back
				x[t][c] = x[t-2][c];
			}
		}
		// Compute the predictive information for block length 1 and check that it
		//  gives us zero bits, since there is no one step correlation
		PredictiveInformationCalculatorDiscrete piCalc = new PredictiveInformationCalculatorDiscrete(2, 1);
		piCalc.initialise();
		piCalc.addObservations(x);
		double pi = piCalc.computeAverageLocalOfObservations();
		assertEquals(0, pi, 0.000000001);
		// Now compute the predictive information for block length 2 and check that it
		//  gives us *two* bits
		piCalc = new PredictiveInformationCalculatorDiscrete(2, 2);
		piCalc.initialise();
		piCalc.addObservations(x);
		pi = piCalc.computeAverageLocalOfObservations();
		assertEquals(2, pi, 0.000000001);
	}
}
