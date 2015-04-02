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

import infodynamics.measures.continuous.MultiInfoCalculator;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.NeighbourNodeData;

/**
 * <p>Computes the differential multi-information of two given multivariate
 *  sets of
 *  observations (implementing {@link MultiInfoCalculator}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see Kraskov et al., below),
 *  <b>algorithm 1</b>.
 *  Most of the functionality is defined by the parent class 
 *  {@link MultiInfoCalculatorKraskov}.</p>
 *
 * <p>Usage is as per the paradigm outlined for {@link MultiInfoCalculator},
 * and expanded on in {@link MultiInfoCalculatorKraskov}.
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
public class MultiInfoCalculatorKraskov1
	extends MultiInfoCalculatorKraskov {
	
	public MultiInfoCalculatorKraskov1() {
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
		
		// Constants:
		double dimensionsMinus1TimesDiGammaN = (double) (dimensions - 1) * digammaN;

		// Count the average number of points within eps_x for each marginal x of each point
		double sumDiGammas = 0;
		double[] sumNMarginals = new double[dimensions];
				
		for (int t = startTimePoint; t < startTimePoint + numTimePoints; t++) {
			// Compute eps for this time step by
			//  finding the kth closest neighbour for point t:
			PriorityQueue<NeighbourNodeData> nnPQ =
					kdTreeJoint.findKNearestNeighbours(k, t, dynCorrExclTime);
			// First element in the PQ is the kth NN,
			//  and epsilon = kthNnData.distance
			NeighbourNodeData kthNnData = nnPQ.poll();
			
			// Count the number of points whose x distance is less
			//  than kthNnData.distance, for each marginal variable x:
			int[] n_marginals = new int[dimensions];
			double thisSumDiGammas = 0;
			for (int d = 0; d < dimensions; d++) {
				n_marginals[d] =
						rangeSearchersInMarginals[d].countPointsStrictlyWithinR(
								t, kthNnData.distance, dynCorrExclTime);
				sumNMarginals[d] += n_marginals[d];
				// And take the digammas:
				thisSumDiGammas += MathsUtils.digamma(n_marginals[d]+1);
			}
			sumDiGammas += thisSumDiGammas;

			if (returnLocals) {
				localMi[t-startTimePoint] = digammaK - thisSumDiGammas +
						dimensionsMinus1TimesDiGammaN;
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
			double[] returnArray = new double[dimensions+1];
			returnArray[0] = sumDiGammas;
			System.arraycopy(sumNMarginals, 0, returnArray, 1, dimensions);
			return returnArray;
		}
	}	
}
