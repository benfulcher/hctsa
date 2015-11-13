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

import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;
import infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov;

/**
 * 
 * = Example 9 - Transfer entropy on continuous data using Kraskov estimators with auto-embedding =
 * 
 * Transfer entropy (TE) calculation on continuous-valued data using the Kraskov-estimator TE calculator,
 * with automatic selection of embedding parameters 
 * 
 * @author Joseph Lizier
 *
 */
public class Example9TeKraskovAutoEmbedding {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		
		// Examine the heart-breath interaction that Schreiber originally looked at:
		ArrayFileReader afr = new ArrayFileReader("../data/SFI-heartRate_breathVol_bloodOx.txt");
		double[][] data = afr.getDouble2DMatrix();
		// Select data points 2350:3550
		data = MatrixUtils.selectRows(data, 2349, 3550-2350+1);

		// Create a Kraskov TE calculator:
		TransferEntropyCalculatorKraskov teCalc = new TransferEntropyCalculatorKraskov();
		
		// Set properties for auto-embedding of both source and destination
		//  using the Ragwitz criteria:
		//  a. Auto-embedding method
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_AUTO_EMBED_METHOD,
				TransferEntropyCalculatorKraskov.AUTO_EMBED_METHOD_RAGWITZ);
		//  b. Search range for embedding dimension (k) and delay (tau)
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_K_SEARCH_MAX, "6");
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_TAU_SEARCH_MAX, "6");
		// Since we're auto-embedding, no need to supply k, l, k_tau, l_tau here:
		teCalc.initialise();
		// Compute TE from breath (column 1) to heart (column 0) 
		teCalc.setObservations(MatrixUtils.selectColumn(data, 1), MatrixUtils.selectColumn(data, 0));
		double teBreathToHeart = teCalc.computeAverageLocalOfObservations();
		// Check the auto-selected parameters and print out the result:
		int optimisedK = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.K_PROP_NAME));
		int optimisedKTau = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.K_TAU_PROP_NAME));
		int optimisedL = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.L_PROP_NAME));
		int optimisedLTau = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.L_TAU_PROP_NAME));
		System.out.printf("TE(breath->heart) was %.3f nats for (heart embedding:) k=%d," +
				"k_tau=%d, (breath embedding:) l=%d,l_tau=%d optimised via Ragwitz criteria\n",
				teBreathToHeart, optimisedK, optimisedKTau, optimisedL, optimisedLTau);

		// Next, embed the destination only using the Ragwitz criteria:
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_AUTO_EMBED_METHOD,
				TransferEntropyCalculatorKraskov.AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY);
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_K_SEARCH_MAX, "6");
		teCalc.setProperty(TransferEntropyCalculatorKraskov.PROP_TAU_SEARCH_MAX, "6");
		// Since we're only auto-embedding the destination, we supply
		//  source embedding here (to overwrite the auto embeddings from above):
		teCalc.setProperty(TransferEntropyCalculatorKraskov.L_PROP_NAME, "1");
		teCalc.setProperty(TransferEntropyCalculatorKraskov.L_TAU_PROP_NAME, "1");
		// Since we're auto-embedding, no need to supply k and k_tau here:
		teCalc.initialise();
		// Compute TE from breath (column 1) to heart (column 0) 
		teCalc.setObservations(MatrixUtils.selectColumn(data, 1), MatrixUtils.selectColumn(data, 0));
		double teBreathToHeartDestEmbedding = teCalc.computeAverageLocalOfObservations();
		// Check the auto-selected parameters and print out the result:
		optimisedK = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.K_PROP_NAME));
		optimisedKTau = Integer.parseInt(teCalc.getProperty(TransferEntropyCalculatorKraskov.K_TAU_PROP_NAME));
		System.out.printf("TE(breath->heart) was %.3f nats for (heart embedding:) k=%d," + 
				"k_tau=%d, optimised via Ragwitz criteria, plus (breath embedding:) l=1,l_tau=1\n",
				teBreathToHeartDestEmbedding, optimisedK, optimisedKTau);

	}
}
