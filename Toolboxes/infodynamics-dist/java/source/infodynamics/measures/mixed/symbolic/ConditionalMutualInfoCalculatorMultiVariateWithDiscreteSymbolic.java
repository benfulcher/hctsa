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

package infodynamics.measures.mixed.symbolic;

import infodynamics.measures.mixed.ConditionalMutualInfoCalculatorMultiVariateWithDiscrete;
import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.measures.discrete.ConditionalMutualInformationCalculatorDiscrete;
import infodynamics.utils.FirstIndexComparatorDouble;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

/**
 * <p>These calculators are <b>EXPERIMENTAL</b> -- not properly tested,
 * and not well documented. The intended calling pattern is similar to
 * {@link ConditionalMutualInfoCalculatorMultiVariate}
 * </p>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSymbolic implements
		ConditionalMutualInfoCalculatorMultiVariateWithDiscrete {

	// The calculator used to do the grunt work
	protected ConditionalMutualInformationCalculatorDiscrete condMiCalc;
	
	protected int dimensions;
	
	protected int[][] permutations;
	// For each permutation index, holds the unique permutation id
	protected int[] permutationIds;
	// For each possible permutation id, holds the permutation index 
	protected int[] idToPermutationIndex;
	
	// Array indices for the 2D array sorted by the first index
	protected static final int VAL_COLUMN = 0;
	protected static final int VAR_NUM_COLUMN = 1;

	public static final String PROP_NORMALISE = "NORMALISE";
	private boolean normalise = true;

	public ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSymbolic() {
		// Nothing to do
	}

	public void initialise(int dimensions, int base, int condBase) throws Exception {
		this.dimensions = dimensions;

		// First work out how many permutations of orderings we could have
		RandomGenerator rg = new RandomGenerator();
		permutations = rg.generateAllDistinctPerturbations(dimensions);
		// Now generate an int signature for each permutation:
		permutationIds = new int[permutations.length];
		for (int r = 0; r < permutations.length; r++) {
			permutationIds[r] = generatePermutationId(permutations[r]);
		}
		// Now we have a list of permutations, each with an ID (which will be dimensions^dimensions)
		// Generate a reverse mapping from permutation identifier to permutation id
		idToPermutationIndex = new int[MathsUtils.power(dimensions, dimensions)];
		// First initialise all mappings to -1 : this will force an Array lookup error if 
		//  an identifier is not calculated correctly (most of the time)
		for (int i = 0; i < idToPermutationIndex.length; i++) {
			idToPermutationIndex[i] = -1;
		}
		for (int idIndex = 0; idIndex < permutationIds.length; idIndex++) {
			idToPermutationIndex[permutationIds[idIndex]] = idIndex;
		}
		
		condMiCalc = ConditionalMutualInformationCalculatorDiscrete.newInstance(permutationIds.length, base, condBase);
		condMiCalc.initialise();
	}

	/**
	 * Generate the unique permutation id for this permutation.
	 * 
	 * @param data
	 * @return
	 */
	private int generatePermutationId(int[] data) {
		int permutationId = 0;
		for (int c = 0; c < dimensions; c++) {
			permutationId *= dimensions;
			permutationId +=  data[c];
		}
		return permutationId;
	}

	/**
	 * Generate the unique permutation id for this permutation.
	 * Convert the floating point variable numbers into ints first
	 * 
	 * @param data
	 * @return
	 */
	private int generatePermutationId(double[] data) {
		int permutationId = 0;
		for (int c = 0; c < dimensions; c++) {
			permutationId *= dimensions;
			permutationId +=  (int) data[c];
		}
		return permutationId;
	}

	public void setProperty(String propertyName, String propertyValue) {
		if (propertyName.equals(PROP_NORMALISE)) {
			normalise = Boolean.parseBoolean(propertyValue);
		}
	}

	public void setObservations(double[][] continuousObservations,
			int[] discreteObservations, int[] conditionalObservations) throws Exception {
		if (normalise) {
			// Normalise the continuous observations first
			continuousObservations = MatrixUtils.normaliseIntoNewArray(continuousObservations);
		}
		// Construct the orderings for the continuous observations
		int[] mappedPermutationIds = new int[continuousObservations.length];
		for (int t = 0; t < continuousObservations.length; t++) {
			// Work out what the order of the continuous variables was here:
			double[][] variablesAndIndices = new double[dimensions][2];
			for (int v = 0; v < dimensions; v++) {
				variablesAndIndices[v][VAL_COLUMN] = continuousObservations[t][v];
				variablesAndIndices[v][VAR_NUM_COLUMN] = v;
			}
			java.util.Arrays.sort(variablesAndIndices, FirstIndexComparatorDouble.getInstance());
			// Now the second column contains the order of values here
			double[] permutation = MatrixUtils.selectColumn(variablesAndIndices, VAR_NUM_COLUMN);
			int permutationId = generatePermutationId(permutation);
			mappedPermutationIds[t] = idToPermutationIndex[permutationId];
		}
		// Now we can set the observations
		condMiCalc.addObservations(mappedPermutationIds, discreteObservations, conditionalObservations);
	}

	public double computeAverageLocalOfObservations() throws Exception {
		return condMiCalc.computeAverageLocalOfObservations();
	}

	public double[] computeLocalUsingPreviousObservations(double[][] contStates,
			int[] discreteStates, int[] conditionedStates) throws Exception {
		throw new Exception("Local method not implemented yet");
	}

	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		return condMiCalc.computeSignificance(numPermutationsToCheck);
	}
	
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings)
		throws Exception {
		return condMiCalc.computeSignificance(newOrderings);
	}

	public void setDebug(boolean debug) {
		// condMiCalc.setDebug(debug);
	}

	public double getLastAverage() {
		return condMiCalc.getLastAverage();
	}

	public int getNumObservations() {
		return condMiCalc.getNumObservations();
	}

}
