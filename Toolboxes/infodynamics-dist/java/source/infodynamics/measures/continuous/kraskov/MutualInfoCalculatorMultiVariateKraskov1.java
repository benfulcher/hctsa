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
 *  <b>algorithm 1</b>.
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
public class MutualInfoCalculatorMultiVariateKraskov1
	extends MutualInfoCalculatorMultiVariateKraskov {
	
	public MutualInfoCalculatorMultiVariateKraskov1() {
		super();
		isAlgorithm1 = true;
	}

	@Override
	protected double[] partialComputeFromObservations(
			int startTimePoint, int numTimePoints, boolean returnLocals) throws Exception {
		
		double startTime = Calendar.getInstance().getTimeInMillis();

		double[] localMi = null;
		if (returnLocals) {
			localMi = new double[numTimePoints];
		}
		
		// Count the average number of points within eps_x and eps_y of each point
		double sumDiGammas = 0;
		double sumNx = 0;
		double sumNy = 0;
				
		for (int t = startTimePoint; t < startTimePoint + numTimePoints; t++) {
			// Compute eps for this time step by
			//  finding the kth closest neighbour for point t:
			PriorityQueue<NeighbourNodeData> nnPQ =
					kdTreeJoint.findKNearestNeighbours(k, t, dynCorrExclTime);
			// First element in the PQ is the kth NN,
			//  and epsilon = kthNnData.distance
			NeighbourNodeData kthNnData = nnPQ.poll();
			
			// Count the number of points whose x distance is less
			//  than eps, and whose y distance is less than
			//  epsilon = kthNnData.distance
			int n_x = nnSearcherSource.countPointsStrictlyWithinR(
							t, kthNnData.distance, dynCorrExclTime);
			int n_y = nnSearcherDest.countPointsStrictlyWithinR(
							t, kthNnData.distance, dynCorrExclTime);
			
			sumNx += n_x;
			sumNy += n_y;
			// And take the digammas:
			double digammaNxPlusOne = MathsUtils.digamma(n_x+1);
			double digammaNyPlusOne = MathsUtils.digamma(n_y+1);
			sumDiGammas += digammaNxPlusOne + digammaNyPlusOne;

			if (returnLocals) {
				// TODO should digammaN be adjusted if we are using 
				//  dynamic correlation exclusion, since we're not sampling from
				//  N after all? Think about this for all Kraskov estimators.
				localMi[t-startTimePoint] = digammaK - digammaNxPlusOne - digammaNyPlusOne + digammaN;
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
