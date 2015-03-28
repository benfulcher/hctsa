/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2015, Joseph T. Lizier
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

package infodynamics.demos.schreiberTransferEntropyExamples;

import infodynamics.utils.ParsedProperties;
import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov;

/**
 * 
 * Used to explore information transfer in the heart rate / breath rate example of Schreiber --
 *  but estimates TE using Kraskov-Stoegbauer-Grassberger estimation.
 *
 * This is intended to be run from the demos/java directory
 *
 * Inputs (via command line arguments)
 * - kHistory - destination embedding length. Defaults to 1.
 * - lHistory - source embedding length. Defaults to 1.
 * - knns - a scalar specifying a single, or X:X format specifying a range, or a comma separated list specifying multiple, value(s) of K nearest neighbours to evaluate TE (Kraskov) with. Defaults to 4.
 * - numSurrogates - a scalar specifying the number of surrogates to evaluate TE from null distribution. Defaults to 0 (i.e. don't evaluate surrogates)
 * Outputs (on standard output), for each knns value:
 * - teHeartToBreath - TE (heart -> breath) for each value of k nearest neighbours
 * - teBreathToHeart - TE (breath -> heart) for each value of k nearest neighbours
 *
 * @author Joseph Lizier
 *
 */
public class HeartBreathRateKraskovRunner {

	/**
	 * @param args a list of kernel widths to evaluate
	 */
	public static void main(String[] args) throws Exception {
		
		// 1. Set up all input parameters:
		int kHistory = 1; // default
		if (args.length >= 1) {
			kHistory = Integer.parseInt(args[0]);
		}
		int lHistory = 1; // default
		if (args.length >= 2) {
			lHistory = Integer.parseInt(args[1]);
		}
		int[] knns;
		if (args.length >= 3) {
			knns = ParsedProperties.parseStringArrayOfInts(args[2]);
		} else {
			// default:
			knns = new int[1];
			knns[0] = 4;
		}
		int numSurrogates = 0; // default
		if (args.length >= 4) {
			numSurrogates = Integer.parseInt(args[3]);
		}
		
		// 2. Load the SFI time-series competition data.
		ArrayFileReader afr = new ArrayFileReader("../data/SFI-heartRate_breathVol_bloodOx.txt");
		double[][] data = afr.getDouble2DMatrix();
		// Select data points 2350:3550
		data = MatrixUtils.selectRows(data, 2349, 3550-2350+1);
		// Separate the data from each column:
		double[] heart = MatrixUtils.selectColumn(data,0); // first column
		double[] chestVol = MatrixUtils.selectColumn(data,1); // second column
		// double[] bloodOx = MatrixUtils.selectColumn(data,2); // Don't need this third column
		int timeSteps = heart.length;

		System.out.printf("TE for heart rate <-> breath rate for kernel estimation with %d samples:\n", timeSteps);

		// Using a KSG transfer entropy estimator
		TransferEntropyCalculatorKraskov teCalc = new TransferEntropyCalculatorKraskov();
		// teCalc.setProperty("NORMALISE", "true"); // this is done by default for this calculator
		
		// Exercise: check if dynamic correlation exclusion (Theiler window) makes a
		//  difference here. Hint: it does change the levels, but not the overall conclusions here:
		// teCalc.setProperty("DYN_CORR_EXCL", "100"); // Dynamic correlation exclusion within 100 samples (in line with what Schreiber did for kernel estimation)
		
		double[] teHeartToBreath = new double[knns.length];
		double[] teBreathToHeart = new double[knns.length];
		
		for (int knnIndex = 0; knnIndex < knns.length; knnIndex++) {
			int knn = knns[knnIndex]; // set next number k of nearest neighbours
			
			// Perform calculation for heart -> breath (lag 1)
			teCalc.initialise(kHistory,1,lHistory,1,1);
			teCalc.setProperty("k", Integer.toString(knn));
			teCalc.setObservations(heart, chestVol);
			teHeartToBreath[knnIndex] = teCalc.computeAverageLocalOfObservations();
			double teHeartToBreathNullMean = 0, teHeartToBreathNullStd = 0;
			if (numSurrogates > 0) {
				EmpiricalMeasurementDistribution teHeartToBreathNullDist = teCalc.computeSignificance(numSurrogates);
				teHeartToBreathNullMean = teHeartToBreathNullDist.getMeanOfDistribution();
				teHeartToBreathNullStd = teHeartToBreathNullDist.getStdOfDistribution();
			}
		
			// Perform calculation for breath -> heart (lag 1)
			teCalc.initialise(kHistory,1,lHistory,1,1);
			teCalc.setProperty("k", Integer.toString(knn));
			teCalc.setObservations(chestVol, heart);
			teBreathToHeart[knnIndex] = teCalc.computeAverageLocalOfObservations();
			double teBreathToHeartNullMean = 0, teBreathToHeartNullStd = 0;
			if (numSurrogates > 0) {
				EmpiricalMeasurementDistribution teBreathToHeartNullDist = teCalc.computeSignificance(numSurrogates);
				teBreathToHeartNullMean = teBreathToHeartNullDist.getMeanOfDistribution();
				teBreathToHeartNullStd = teBreathToHeartNullDist.getStdOfDistribution();
			}
			
			System.out.printf("TE(k=%d,l=%d,knn=%d): h->b = %.3f", kHistory, lHistory, knn, teHeartToBreath[knnIndex]);
			if (numSurrogates > 0) {
				System.out.printf(" (null = %.3f +/- %.3f)", teHeartToBreathNullMean, teHeartToBreathNullStd);
			}
			System.out.printf(", b->h = %.3f nats", teBreathToHeart[knnIndex]);
			if (numSurrogates > 0) {
				System.out.printf("(null = %.3f +/- %.3f)\n", teBreathToHeartNullMean, teBreathToHeartNullStd);
			} else {
				System.out.println();
			}
		}		
	}
}
