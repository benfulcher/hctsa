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

package infodynamics.measures.continuous.symbolic;

import infodynamics.measures.continuous.TransferEntropyCalculator;
import infodynamics.measures.continuous.TransferEntropyCommon;
import infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete;
import infodynamics.utils.FirstIndexComparatorDouble;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.RandomGenerator;

import java.util.Iterator;
import java.util.Vector;

/**
 * <p>Computes the differential transfer entropy (TE) between two univariate
 *  <code>double[]</code> time-series of observations
 *  using symbolic TE estimation.
 *  For details on symbolic TE estimation, see Staniek and Lehnertz (below).</p>
 *  
 * <p>TE was defined by Schreiber (below).
 *  This estimator is realised here by extending
 *  {@link TransferEntropyCommon}.</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link TransferEntropyCalculator},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #TransferEntropyCalculatorSymbolic()}.</li>
 *  <li>A restricted set of properties are available, see {@link #setProperty(String, String)};</li>
 *  </ul>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. Schreiber, <a href="http://dx.doi.org/10.1103/PhysRevLett.85.461">
 * "Measuring information transfer"</a>,
 *  Physical Review Letters 85 (2) pp.461-464, 2000.</li>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href="http://dx.doi.org/10.1103/PhysRevE.77.026110">
 *  "Local information transfer as a spatiotemporal filter for complex systems"</a>
 *  Physical Review E 77, 026110, 2008.</li>
 *  <li>M. Staniek, K. Lehnertz,
 *  <a href="http://dx.doi.org/10.1103/physrevlett.100.158101">
 *  "Symbolic Transfer Entropy"</a>,
 *  Physical Review Letters 100, 158101 (2008)</li>
 * </ul>
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 *
 */
public class TransferEntropyCalculatorSymbolic
		extends TransferEntropyCommon
		implements TransferEntropyCalculator {

	/**
	 * The discrete calculator used to do the grunt work
	 */
	protected TransferEntropyCalculatorDiscrete teCalc;

	protected int maxEmbeddingLength; // = Max(l,k)
	// l will default to the value of k unless l gets explicitly set
	protected int l = 1;
	protected boolean lWasSet = false;
	
	private int[][] destPermutations;
	// For each permutation index, holds the unique permutation id
	private int[] destPermutationIds;
	// For each possible permutation id, holds the permutation index 
	private int[] idToDestPermutationIndex;
	
	private int[][] sourcePermutations;
	// For each permutation index, holds the unique permutation id
	private int[] sourcePermutationIds;
	// For each possible permutation id, holds the permutation index 
	private int[] idToSourcePermutationIndex;

	// Storage for the computation symbols
	protected Vector<int[]> destSymbolsVector;
	protected Vector<int[]> sourceSymbolsVector;
	
	// Array indices for the 2D array sorted by the first index
	protected static final int VAL_COLUMN = 0;
	protected static final int VAR_NUM_COLUMN = 1;

	/**
	 * Construct an instance
	 */
	public TransferEntropyCalculatorSymbolic() {
		// Nothing to do
	}

	@Override
	public void initialise(int k) throws Exception {
		super.initialise(k); // calls initialise();
	}

	@Override
	public void initialise() throws Exception {
		if (!lWasSet) {
			l = k;
		}
		
		// The discrete TE calculator will be run with a base of the number of
		//  permutations of the maximum embedding length (out of l,k), and history 1
		//  since the symbols incorporate the past k states.
		maxEmbeddingLength = Math.max(k, l);
		
		// First work out how many permutations of orderings we could have
		RandomGenerator rg = new RandomGenerator();
		destPermutations = rg.generateAllDistinctPerturbations(k);
		// Now generate an int signature for each permutation:
		destPermutationIds = new int[destPermutations.length];
		for (int r = 0; r < destPermutations.length; r++) {
			destPermutationIds[r] = generatePermutationId(destPermutations[r], k);
		}
		// Now we have a list of permutations, each with an ID (which will be up to k^k)
		// Generate a reverse mapping from permutation identifier to permutation id
		idToDestPermutationIndex = new int[MathsUtils.power(k, k)];
		// First initialise all mappings to -1 : this will force an Array lookup error if 
		//  an identifier is not calculated correctly (most of the time)
		for (int i = 0; i < idToDestPermutationIndex.length; i++) {
			idToDestPermutationIndex[i] = -1;
		}
		for (int idIndex = 0; idIndex < destPermutationIds.length; idIndex++) {
			idToDestPermutationIndex[destPermutationIds[idIndex]] = idIndex;
			/* System.out.print("Permutation "+ idIndex + ": ");
			for (int j = 0; j < permutations[idIndex].length; j++) {
				System.out.print(permutations[idIndex][j] + " ");
			}
			System.out.println(" -> id " + permutationIds[idIndex]);
			*/
		}

		// Now generate the permutations for the source
		sourcePermutations = rg.generateAllDistinctPerturbations(l);
		// Now generate an int signature for each permutation:
		sourcePermutationIds = new int[sourcePermutations.length];
		for (int r = 0; r < sourcePermutations.length; r++) {
			sourcePermutationIds[r] = generatePermutationId(sourcePermutations[r], l);
		}
		// Now we have a list of permutations, each with an ID (which will be up to l^l)
		// Generate a reverse mapping from permutation identifier to permutation id
		idToSourcePermutationIndex = new int[MathsUtils.power(l, l)];
		// First initialise all mappings to -1 : this will force an Array lookup error if 
		//  an identifier is not calculated correctly (most of the time)
		for (int i = 0; i < idToSourcePermutationIndex.length; i++) {
			idToSourcePermutationIndex[i] = -1;
		}
		for (int idIndex = 0; idIndex < sourcePermutationIds.length; idIndex++) {
			idToSourcePermutationIndex[sourcePermutationIds[idIndex]] = idIndex;
			/* System.out.print("Permutation "+ idIndex + ": ");
			for (int j = 0; j < permutations[idIndex].length; j++) {
				System.out.print(permutations[idIndex][j] + " ");
			}
			System.out.println(" -> id " + permutationIds[idIndex]);
			*/
		}

		// The discrete calculator only uses a history of 1 here - k is built into
		//  the permutations
		int base = Math.max(destPermutationIds.length, sourcePermutationIds.length);
		teCalc = TransferEntropyCalculatorDiscrete.newInstance(base, 1);
		teCalc.initialise();
	}

	/**
	 * <p>Set properties for the kernel TE calculator.
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 * 
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 * 		<li>{@link #L_PROP_NAME} -- embedding length for the source past history vector
	 * 			(default value 1)</li>
	 *  	<li>any valid properties for {@link TransferEntropyCommon#setProperty(String, String)},
	 *  		noting that we are explicitly overriding it to allow {@link #L_PROP_NAME} here.</li>
	 * </ul>
	 * </p>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property
	 * @throws Exception for invalid property values
	 */
	@Override
	public void setProperty(String propertyName, String propertyValue) throws Exception {
		boolean propertySet = true;
		if (propertyName.equalsIgnoreCase(L_PROP_NAME)) {
			l = Integer.parseInt(propertyValue);
			lWasSet = true;
		} else {
			// No property was set
			propertySet = false;
			// try the superclass:
			super.setProperty(propertyName, propertyValue);
		}
		if (debug && propertySet) {
			System.out.println("Set property " + propertyName +
					" to " + propertyValue);
		}
	}

	@Override
	public void finaliseAddObservations() {
		// First work out the size to allocate the joint vectors, and do the allocation:
		totalObservations = 0;
		for (double[] destination : vectorOfDestinationObservations) {
			totalObservations += destination.length - maxEmbeddingLength;
		}

		destSymbolsVector = new Vector<int[]>();
		sourceSymbolsVector = new Vector<int[]>();
		
		// Construct the symbols from the given observations
		//int startObservation = 0;
		Iterator<double[]> iterator = vectorOfDestinationObservations.iterator();
		for (double[] source : vectorOfSourceObservations) {
			double[] destination = iterator.next();
			// compute the embedding vectors for the destination and source
			double[][] currentDestPastVectors = null;
			double[][] sourceVectors = null;
			try {
				// we don't use the pre-defined makeDestPastVectors method because it skips off 
				//  the last embedding vector
				currentDestPastVectors = MatrixUtils.makeDelayEmbeddingVector(destination, k, k-1, destination.length - k + 1);
				sourceVectors = MatrixUtils.makeDelayEmbeddingVector(source, l, l-1, source.length - l + 1);
			} catch (Exception e) {
				// The parameters for the above call should be fine, so we don't expect to
				//  throw an Exception here - embed in a RuntimeException if it occurs 
				throw new RuntimeException(e);
			}
			
			// Now compute the permutation values for the dest
			double[][] destVariablesAndIndices = new double[k][2];
			int[] destSymbols = new int[currentDestPastVectors.length];
			for (int t = 0; t < currentDestPastVectors.length; t++) {
				// Work out what the order of the embedded variables was here:
				for (int v = 0; v < k; v++) {
					destVariablesAndIndices[v][VAL_COLUMN] = currentDestPastVectors[t][v];
					destVariablesAndIndices[v][VAR_NUM_COLUMN] = v;
				}
				java.util.Arrays.sort(destVariablesAndIndices, FirstIndexComparatorDouble.getInstance());
				// Now the second column contains the order of values here
				double[] permutation = MatrixUtils.selectColumn(destVariablesAndIndices, VAR_NUM_COLUMN);
				int permutationId = generatePermutationId(permutation, k);
				destSymbols[t] = idToDestPermutationIndex[permutationId];
			}

			// Compute the permutation values for the source
			double[][] sourceVariablesAndIndices = new double[l][2];
			int[] sourceSymbols = new int[sourceVectors.length];
			for (int t = 0; t < sourceVectors.length; t++) {
				// Work out what the order of the embedded variables was here:
				for (int v = 0; v < l; v++) {
					sourceVariablesAndIndices[v][VAL_COLUMN] = sourceVectors[t][v];
					sourceVariablesAndIndices[v][VAR_NUM_COLUMN] = v;
				}
				java.util.Arrays.sort(sourceVariablesAndIndices, FirstIndexComparatorDouble.getInstance());
				// Now the second column contains the order of values here
				double[] permutation = MatrixUtils.selectColumn(sourceVariablesAndIndices, VAR_NUM_COLUMN);
				int permutationId = generatePermutationId(permutation, l);
				sourceSymbols[t] = idToSourcePermutationIndex[permutationId];
			}
			
			if (k > l) {
				// l is smaller, so we will have more source embeddings than destination ones.
				// Get rid of the first k - l source embeddings
				sourceSymbols = MatrixUtils.select(sourceSymbols, k - l,
						sourceSymbols.length - (k - l));
			} else if (l > k) {
				// k is smaller, so we will have more dest embeddings than source ones.
				// Get rid of the first l - k source embeddings
				destSymbols = MatrixUtils.select(destSymbols, l - k,
						destSymbols.length - (l - k));
			}
			
			// Add these observations in, and keep them locally
			destSymbolsVector.add(destSymbols);
			sourceSymbolsVector.add(sourceSymbols);
			teCalc.addObservations(sourceSymbols, destSymbols);
			
			if (destination.length - maxEmbeddingLength != destSymbols.length - 1) {
				throw new RuntimeException(
					String.format("Number of observations %d doesn't match what's expected %d",
						destSymbols.length - 1, destination.length - maxEmbeddingLength));
			}
			
			// I don't remember why we were tracking this:
			// (probably before using vectors of observations)
			//startObservation += destSymbols.length - 1;
		}
		
	}

	/**
	 * Generate the unique permutation id for this permutation.
	 * 
	 * @param data
	 * @param base the number of elements being permuted
	 * @return
	 */
	private int generatePermutationId(int[] data, int base) {
		int permutationId = 0;
		for (int c = 0; c < data.length; c++) {
			permutationId *= base;
			permutationId +=  data[c];
		}
		return permutationId;
	}

	/**
	 * Generate the unique permutation id for this ordering of variable ids.
	 * Convert the floating point variable ids into ints first
	 * 
	 * @param ids
	 * @param base the number of elements being permuted
	 * @return
	 */
	private int generatePermutationId(double[] ids, int base) {
		int permutationId = 0;
		for (int c = 0; c < ids.length; c++) {
			permutationId *= base;
			permutationId +=  (int) ids[c];
		}
		return permutationId;
	}

	@Override
	public double computeAverageLocalOfObservations() throws Exception {
		lastAverage = teCalc.computeAverageLocalOfObservations();
		return lastAverage;
	}

	@Override
	public double[] computeLocalOfPreviousObservations() throws Exception {
		double[] locals = new double[totalObservations];
		int currentIndexInLocals = 0;
		Iterator<int[]> iterator = destSymbolsVector.iterator();
		for (int[] sourceSymbols : sourceSymbolsVector) {
			int[] destSymbols = iterator.next();
			double[] theseLocals = teCalc.computeLocalFromPreviousObservations(sourceSymbols, destSymbols);
			System.arraycopy(theseLocals, 1, locals, currentIndexInLocals, theseLocals.length - 1);
			currentIndexInLocals += theseLocals.length - 1;
		}
		lastAverage = MatrixUtils.mean(locals);
		return locals;
	}

	/**
	 * Not implemented yet
	 */
	@Override
	public double[] computeLocalUsingPreviousObservations(
			double[] newSourceObservations, double[] newDestObservations)
			throws Exception {
		throw new UnsupportedOperationException("Not implemented yet");
	}

	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int numPermutationsToCheck) throws Exception {
		return teCalc.computeSignificance(numPermutationsToCheck);
	}

	/**
	 * Not implemented yet
	 */
	public EmpiricalMeasurementDistribution computeSignificance(int[][] newOrderings)
			throws Exception {
		throw new UnsupportedOperationException("Not implemented yet");
		/* This method used to cheat:
		System.out.println("TESymbolic.computeSignificance(): Not using the new orderings supplied");
		return teCalc.computeSignificance(newOrderings.length);
		*/
	}
}
