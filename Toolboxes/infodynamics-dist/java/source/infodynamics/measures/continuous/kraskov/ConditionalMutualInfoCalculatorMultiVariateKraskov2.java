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

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.NeighbourNodeData;

/**
 * <p>Computes the differential conditional mutual information of two multivariate
 *  <code>double[][]</code> sets of observations, conditioned on another
 *  (implementing {@link ConditionalMutualInfoCalculatorMultiVariate}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see references below)
 *  <b>algorithm 2</b>.
 *  Most of the functionality is defined by the parent class 
 *  {@link ConditionalMutualInfoCalculatorMultiVariateKraskov}.</p>
 *
 * <p>Crucially, the calculation is performed by examining
 * neighbours in the full joint space
 * rather than two MI calculators.
 * This is roughly as specified by Frenzel and Pompe (who only
 * specified this for algorithm 1), but mathematically
 * adapted to algorithm 2 by Wibral et al. (see below).</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link ConditionalMutualInfoCalculatorMultiVariate},
 * and expanded on in {@link ConditionalMutualInfoCalculatorMultiVariateKraskov}.
 * </p>
 *  
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>M. Wibral, R. Vicente, and M. Lindner, <a href="http://dx.doi.org/10.1007/978-3-642-54474-3_1">
 *  "Transfer Entropy in Neuroscience"</a>,
 *  in "Directed Information Measures in Neuroscience",
 *  Understanding Complex Systems series, edited by M. Wibral, R. Vicente, and J. T. Lizier
 *  (Springer, Berlin/Heidelberg, 2014) pp. 3--36.</li>
 * 	<li>Frenzel and Pompe, <a href="http://dx.doi.org/10.1103/physrevlett.99.204101">
 * 	"Partial Mutual Information for Coupling Analysis of Multivariate Time Series"</a>,
 * 	Physical Review Letters, <b>99</b>, p. 204101+ (2007).</li>
 * 	<li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @author Ipek Özdemir
 */
public class ConditionalMutualInfoCalculatorMultiVariateKraskov2
	extends ConditionalMutualInfoCalculatorMultiVariateKraskov {

	@Override
	protected double[] partialComputeFromObservations(int startTimePoint,
			int numTimePoints, boolean returnLocals) throws Exception {
		
		double startTime = Calendar.getInstance().getTimeInMillis();

		double[] localCondMi = null;
		if (returnLocals) {
			localCondMi = new double[numTimePoints];
		}
		
		// Constants:
		double twoOnK = 2.0 / (double) k;
		
		// Count the average number of points within eps_xz, eps_yz and eps_z
		double sumDiGammas = 0;
		double sumNxz = 0;
		double sumNyz = 0;
		double sumNz = 0;
		double sumInverseCountInJointYZ = 0;
		double sumInverseCountInJointXZ = 0;

		// Arrays used for fast searching on conditionals with a marginal:
		boolean[] isWithinRForConditionals = new boolean[totalObservations];
		int[] indicesWithinRForConditionals = new int[totalObservations+1];
		
		for (int t = startTimePoint; t < startTimePoint + numTimePoints; t++) {
			// Compute eps_x and eps_y and eps_z for this time step by
			//  finding the kth closest neighbours for point t:
			PriorityQueue<NeighbourNodeData> nnPQ =
					kdTreeJoint.findKNearestNeighbours(k, t, dynCorrExclTime);
			
			// Find eps_{x,y,z} as the maximum x, y and z norms amongst this set:
			double eps_x = 0.0;
			double eps_y = 0.0;
			double eps_z = 0.0;
			for (int j = 0; j < k; j++) {
				// Take the furthest remaining of the nearest neighbours from the PQ:
				NeighbourNodeData nnData = nnPQ.poll();
				if (nnData.norms[0] > eps_x) {
					eps_x = nnData.norms[0];
				}
				if (nnData.norms[1] > eps_y) {
					eps_y = nnData.norms[1];
				}
				if (nnData.norms[2] > eps_z) {
					eps_z = nnData.norms[2];
				}
			}

			// Count the number of points whose z distance is less
			//  than or equal to eps_z, and whose x and z distances are less
			//  than or equal to eps_z and eps_x, and whose y and z distance are less
			//  than or equal to eps_z and eps_y:

			/* Option A -- straightforward way using each k-d tree separately:
			int n_xz = nnSearcherVar1.countPointsWithinOrOnRs(
					t, new double[] {eps_x, eps_z}, dynCorrExclTime);
			int n_yz = nnSearcherVar2.countPointsWithinOrOnRs(
					t, new double[] {eps_y, eps_z}, dynCorrExclTime);
			int n_z = nnSearcherConditional.countPointsWithinOrOnR(
					t, eps_z, dynCorrExclTime);
			*/
			
			// Option C --
			// Identify the points satisfying the conditional criteria, then use
			//  the knowledge of which points made this cut to speed up the searching
			//  in the conditional-marginal spaces:
			// 1. Identify the n_z points within the conditional boundaries:
			nnSearcherConditional.findPointsWithinR(t, eps_z, dynCorrExclTime,
					true, isWithinRForConditionals, indicesWithinRForConditionals);
			// 2. Then compute n_xz and n_yz harnessing our knowledge of
			//  which points qualified for the conditional already:
			// Don't need to supply dynCorrExclTime in the following, because only 
			//  points outside of it have been included in isWithinRForConditionals
			int n_xz;
			if (dimensionsVar1 > 1) {
				// Check only the x variable against eps_x, use existing results for z
				n_xz = kdTreeVar1Conditional.countPointsWithinRs(t, 
						new double[] {eps_x, eps_z},
						true, 1, isWithinRForConditionals);
			} else { // Generally faster to search only the marginal space if it is univariate 
				n_xz = uniNNSearcherVar1.countPointsWithinR(t, eps_x,
						true, isWithinRForConditionals);
			}
			int n_yz;
			if (dimensionsVar2 > 1) {
				// Check only the y variable against eps_y, use existing results for z
				n_yz = kdTreeVar2Conditional.countPointsWithinRs(t,
						new double[] {eps_y, eps_z},
						true, 1, isWithinRForConditionals);
			} else { // Generally faster to search only the marginal space if it is univariate
				n_yz = uniNNSearcherVar2.countPointsWithinR(t, eps_y,
						true, isWithinRForConditionals);
			}
			// 3. Finally, reset our boolean array for its next use while we count n_z:
			int n_z;
			for (n_z = 0; indicesWithinRForConditionals[n_z] != -1; n_z++) {
				isWithinRForConditionals[indicesWithinRForConditionals[n_z]] = false;
			}
			// end option C

			
			sumNxz += n_xz;
			sumNyz += n_yz;
			sumNz += n_z;
			// And take the digammas:
			double digammaNxz = MathsUtils.digamma(n_xz);
			double digammaNyz = MathsUtils.digamma(n_yz);
			double digammaNz = MathsUtils.digamma(n_z);
			double invN_xz = 1.0/(double) n_xz;
			double invN_yz = 1.0/(double) n_yz;
			sumInverseCountInJointXZ += invN_xz;
			sumInverseCountInJointYZ += invN_yz;
			double contributionDigammas = digammaNz - digammaNxz - digammaNyz;
			sumDiGammas += contributionDigammas;

			if (returnLocals) {
				localCondMi[t-startTimePoint] = digammaK - twoOnK + contributionDigammas + invN_xz + invN_yz;
				if (debug) {
					// Only tracking this for debugging purposes: 
					System.out.printf("t=%d, n_xz=%d, n_yz=%d, n_z=%d, 1/n_yz=%.3f, 1/n_xz=%.3f, local=%.4f\n",
						t, n_xz, n_yz, n_z, invN_yz, invN_xz, localCondMi[t-startTimePoint]);
				}
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
			return localCondMi;
		} else {
			double[] returnValues = new double[6];
			returnValues[0] = sumDiGammas;
			returnValues[1] = sumNxz;
			returnValues[2] = sumNyz;
			returnValues[3] = sumNz;
			returnValues[4] = sumInverseCountInJointXZ;
			returnValues[5] = sumInverseCountInJointYZ;
			return returnValues;
		}		
	}
}
