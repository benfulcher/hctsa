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

package infodynamics.networkinference.interregional;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;

/**
 * Extends MeasurementDistribution for computations over 
 *  large sets where we generate measures over permutations of subsets.
 * 
 * The member distribution refers to averages over subsets for each permutation.
 * The new member distributionForSubsets refers to averages over permutations
 *   for each subset.
 * 
 * @author Joseph Lizier
 *
 */
public class MeasurementDistributionPermutationsOverSubsets extends EmpiricalMeasurementDistribution {

	// The true measurements for each subset s
	double[] actualValues;
	
	// Distribution over [subsets] then [permutations]
	double[][] distributionOverSubsetsAndPermutations;
	
	// The average over permutations for each subset s
	double[] avDistributionForSubsets;
	
	/**
	 * 
	 */
	public MeasurementDistributionPermutationsOverSubsets(int size) {
		super(size);
	}

	public MeasurementDistributionPermutationsOverSubsets
		(double[][] theDistribution, double[] actualValues) {

		super(theDistribution[0].length);

		this.actualValues = actualValues;
		this.actualValue = MatrixUtils.mean(actualValues);
		distributionOverSubsetsAndPermutations = theDistribution;
		
		int subsets = theDistribution.length;
		int permutations = theDistribution[0].length;
		// distribution, which is av over permutations, was created in the
		//  super constructor.
		avDistributionForSubsets = new double[subsets];
		
		// distribution will hold the averages for each permutation i
		int avValuesFromDistributionGreaterThanActualAvs = 0;
		for (int i = 0; i < permutations; i++) {
			distribution[i] = MatrixUtils.mean(theDistribution, i);
			if (distribution[i] >= actualValue) {
				avValuesFromDistributionGreaterThanActualAvs++;
			}
		}
		pValue = 
			(double) avValuesFromDistributionGreaterThanActualAvs / (double) permutations;
		
		for (int s = 0; s < subsets; s++) {
			avDistributionForSubsets[s] = MatrixUtils.mean(theDistribution[s]);
		}
	}
}
