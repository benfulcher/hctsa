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

package infodynamics.measures.continuous.kraskov;

import java.util.Calendar;
import java.util.PriorityQueue;

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.NeighbourNodeData;

/**
 * <p>Computes the differential mutual information of two given multivariate sets of
 *  observations (implementing {@link MutualInfoCalculatorMultiVariate}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see Kraskov et al., below),
 *  <b>algorithm 2</b>.
 *  Most of the functionality is defined by the parent class 
 *  {@link MutualInfoCalculatorMultiVariateKraskov}.</p>
 *
 * <p>Usage is as per the paradigm outlined for {@link MutualInfoCalculatorMultiVariate},
 * and expanded on in {@link MutualInfoCalculatorMultiVariateKraskov}.
 * </p>
 *  
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @author Ipek Özdemir
 */
public class MutualInfoCalculatorMultiVariateKraskov2
	extends MutualInfoCalculatorMultiVariateKraskov {
	
	public MutualInfoCalculatorMultiVariateKraskov2() {
		super();
		isAlgorithm1 = false;
	}

	protected double[] partialComputeFromObservations(
			int startTimePoint, int numTimePoints, boolean returnLocals) throws Exception {
		
		double startTime = Calendar.getInstance().getTimeInMillis();

		double[] localMi = null;
		if (returnLocals) {
			localMi = new double[numTimePoints];
		}
		
		// Constants:
		double invK = 1.0 / (double)k;
		
		// Count the average number of points within eps_x and eps_y of each point
		double sumDiGammas = 0;
		double sumNx = 0;
		double sumNy = 0;
				
		for (int t = startTimePoint; t < startTimePoint + numTimePoints; t++) {
			// Compute eps_x and eps_y for this time step by
			//  finding the kth closest neighbours for point t:
			PriorityQueue<NeighbourNodeData> nnPQ =
					kdTreeJoint.findKNearestNeighbours(k, t, dynCorrExclTime);

			// Find eps_{x,y} as the maximum x and y norms amongst this set:
			double eps_x = 0.0;
			double eps_y = 0.0;
			for (int j = 0; j < k; j++) {
				// Take the furthest remaining of the nearest neighbours from the PQ:
				NeighbourNodeData nnData = nnPQ.poll();
				if (nnData.norms[0] > eps_x) {
					eps_x = nnData.norms[0];
				}
				if (nnData.norms[1] > eps_y) {
					eps_y = nnData.norms[1];
				}
			}
			
			// Count the number of points whose x distance is less
			//  than or equal to eps_x, and whose y distance is less
			//  than or equal to eps_y
			int n_x = nnSearcherSource.countPointsWithinOrOnR(
					t, eps_x, dynCorrExclTime);
			int n_y = nnSearcherDest.countPointsWithinOrOnR(
					t, eps_y, dynCorrExclTime);

			sumNx += n_x;
			sumNy += n_y;
			// And take the digammas:
			double digammaNx = MathsUtils.digamma(n_x);
			double digammaNy = MathsUtils.digamma(n_y);
			sumDiGammas += digammaNx + digammaNy;

			if (returnLocals) {
				localMi[t-startTimePoint] = digammaK - invK - digammaNx - digammaNy + digammaN;
			}
		}

		if (debug) {
			Calendar rightNow2 = Calendar.getInstance();
			long endTime = rightNow2.getTimeInMillis();
			System.out.println("Subset " + startTimePoint + ":" +
					(startTimePoint + numTimePoints) + " Calculation time: " +
					((endTime - startTime)/1000.0) + " sec" );
		}
		
		// Select what to return:
		if (returnLocals) {
			return localMi;
		} else {
			return new double[] {sumDiGammas, sumNx, sumNy};
		}
	}
}
