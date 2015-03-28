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

import infodynamics.utils.MatrixUtils;

/**
 * Implements <b>separable information</b> (see Lizier et al, 2010, below),
 *  by using separate Transfer Entropy and Active information storage
 *  calculators.
 *  
 * Javadocs are preliminary here (TODO), please see 
 * {@link SeparableInfoCalculatorDiscrete} for user-level documentation
 * regarding the methods.
 * 
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href=http://dx.doi.org/10.1063/1.3486801">
 *  "Information modification and particle collisions in distributed computation"</a>
 *  Chaos 20, 3, 037109 (2010).</li>
 * </ul>
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class SeparableInfoCalculatorDiscreteByAddition extends SeparableInfoCalculatorDiscrete {

	ActiveInformationCalculatorDiscrete aiCalc;
	TransferEntropyCalculatorDiscrete[] ateCalcs;
	
	private boolean localStatsValid = true;
	
	public SeparableInfoCalculatorDiscreteByAddition(int base, int history,
			int numInfoContributors) {
		// Create super class without creating any storage
		//  for the observations
		super(base, history, numInfoContributors, true);
		
		// Create the calculators
		aiCalc = ActiveInformationCalculatorDiscrete.newInstance(base, history);
	}

	/**
	 * Explicitly create the transfer entropy calculators
	 *  when required (this can save much memory in certain circumstances,
	 *  e.g. calling computeAverageLocal methods).
	 *
	 */
	private void createAppTeCalculators() {
		ateCalcs = new TransferEntropyCalculatorDiscrete[numSources];
		for (int i = 0; i < numSources; i++) {
			ateCalcs[i] = TransferEntropyCalculatorDiscrete.newInstance(base, k);
			ateCalcs[i].setPeriodicBoundaryConditions(periodicBoundaryConditions);
		}
	}
	
	@Override
	public void initialise() {
		super.initialise();
		aiCalc.initialise();
		if (ateCalcs == null) {
			createAppTeCalculators();
		}
		for (int i = 0; i < numSources; i++) {
			ateCalcs[i].initialise();
		}		
	}

	@Override
	public void addObservations(int[][] states, int destCol, int[] sourcesAbsolute) {
		aiCalc.addObservations(states);
		int[] cleanedSourcesAbsolute = cleanAbsoluteSources(sourcesAbsolute, destCol);
		for (int i = 0; i < numSources; i++) {
			ateCalcs[i].addObservations(states, cleanedSourcesAbsolute[i], destCol);
		}
	}

	@Override
	public void addObservations(int[][] states, int[] offsetOfDestFromSources) {
		aiCalc.addObservations(states);
		int[] cleanedOffsets = cleanOffsetOfDestFromSources(offsetOfDestFromSources);
		for (int i = 0; i < numSources; i++) {
			ateCalcs[i].addObservations(states, cleanedOffsets[i]);
		}
	}

	@Override
	public void addObservations(int[][][] states, int destAgentRow, int destAgentColumn, int[][] sourcesAbsolute) {
		aiCalc.addObservations(states, destAgentRow, destAgentColumn);
		int[][] cleanedSourcesAbsolute = cleanAbsoluteSources(sourcesAbsolute, destAgentRow, destAgentColumn);
		for (int i = 0; i < numSources; i++) {
			ateCalcs[i].addObservations(states, cleanedSourcesAbsolute[i][ROW_INDEX], cleanedSourcesAbsolute[i][COLUMN_INDEX],
					destAgentRow, destAgentColumn);
		}
	}

	@Override
	public void addObservations(int[][][] states, int[][] offsetOfDestFromSources) {
		aiCalc.addObservations(states);
		int[][] cleanedOffsets = cleanOffsetOfDestFromSources(offsetOfDestFromSources);
		for (int i = 0; i < numSources; i++) {
			ateCalcs[i].addObservations(states, cleanedOffsets[i][ROW_INDEX],
					cleanedOffsets[i][COLUMN_INDEX]);
		}
	}

	@Override
	public synchronized double computeAverageLocalOfObservations() {
		average = aiCalc.computeAverageLocalOfObservations();
		for (int i = 0; i < numSources; i++) {
			average += ateCalcs[i].computeAverageLocalOfObservations();
		}
		localStatsValid = false;
		return average;
	}

	@Override
	public double[] computeLocalFromPreviousObservations(int[][] states, int destCol, int[] sourcesAbsolute) {
		double[] local = aiCalc.computeLocalFromPreviousObservations(states, destCol);
		double[][] localComponents = null;
		if (computeMultiInfoCoherence) {
			localComponents = new double[numSources + 1][];
			// Keep a link to the local active info
			localComponents[0] = local;
		}
		int[] cleanedSourcesAbsolute = cleanAbsoluteSources(sourcesAbsolute, destCol);
		double[] temp;
		for (int i = 0; i < numSources; i++) {
			temp = ateCalcs[i].computeLocalFromPreviousObservations(states, cleanedSourcesAbsolute[i], destCol);
			if (computeMultiInfoCoherence) {
				// Keep a link to this apparent transfer
				localComponents[1 + numSources] = temp;
			}
			try {
				local = MatrixUtils.add(local, temp);
			} catch (Exception e) {
				// Exception only thrown where arrays were not
				//  of the same length - should not happen here.
				return null;
			}
		}
		// Set statistics and if not periodic boundary conditions, 
		//  clean up points which don't get all information contributors
		setStatistics(local, localComponents);
		localStatsValid = true;
		return local;
	}

	@Override
	public double[][] computeLocalFromPreviousObservations(int[][] states, int[] offsetOfDestFromSources) {
		double[][] local = aiCalc.computeLocalFromPreviousObservations(states);
		double[][][] localComponents = null;
		if (computeMultiInfoCoherence) {
			localComponents = new double[numSources + 1][][];
			// Keep a link to the local active info
			localComponents[0] = local;
		}
		int[] cleanedOffsets = cleanOffsetOfDestFromSources(offsetOfDestFromSources);
		double[][] temp;
		for (int i = 0; i < numSources; i++) {
			temp = ateCalcs[i].computeLocalFromPreviousObservations(
					states, cleanedOffsets[i]);
			if (computeMultiInfoCoherence) {
				// Keep a link to this apparent transfer
				localComponents[1 + numSources] = temp;
			}
			try {
				local = MatrixUtils.add(local, temp);
			} catch (Exception e) {
				// Exception only thrown where arrays were not
				//  of the same length - should not happen here.
				return null;
			}
		}
		// Set statistics and if not periodic boundary conditions, 
		//  clean up points which don't get all information contributors
		setStatistics(local, cleanedOffsets, localComponents);
		localStatsValid = true;
		return local;
	}

	@Override
	public double[] computeLocalFromPreviousObservations(int[][][] states,
			int destAgentRow, int destAgentColumn, int[][] sourcesAbsolute) {
		double[] local = aiCalc.computeLocalFromPreviousObservations(states, destAgentRow, destAgentColumn);
		double[][] localComponents = null;
		if (computeMultiInfoCoherence) {
			localComponents = new double[numSources + 1][];
			// Keep a link to the local active info
			localComponents[0] = local;
		}
		int[][] cleanedSourcesAbsolute = cleanAbsoluteSources(sourcesAbsolute, destAgentRow, destAgentColumn);
		double[] temp;
		for (int i = 0; i < numSources; i++) {
			temp = ateCalcs[i].computeLocalFromPreviousObservations(
					states, cleanedSourcesAbsolute[i][ROW_INDEX], cleanedSourcesAbsolute[i][COLUMN_INDEX],
					destAgentRow,
					destAgentColumn);
			if (computeMultiInfoCoherence) {
				// Keep a link to this apparent transfer
				localComponents[1 + numSources] = temp;
			}
			try {
				local = MatrixUtils.add(local, temp);
			} catch (Exception e) {
				// Exception only thrown where arrays were not
				//  of the same length - should not happen here.
				return null;
			}
		}
		// Set statistics
		setStatistics(local, localComponents);
		localStatsValid = true;
		return local;
	}

	@Override
	public double[][][] computeLocalFromPreviousObservations(int[][][] states, int[][] offsetOfDestFromSources) {
		double[][][] local = aiCalc.computeLocalFromPreviousObservations(states);
		double[][][][] localComponents = null;
		if (computeMultiInfoCoherence) {
			localComponents = new double[numSources + 1][][][];
			// Keep a link to the local active info
			localComponents[0] = local;
		}
		int[][] cleanedOffsets = cleanOffsetOfDestFromSources(offsetOfDestFromSources);
		double[][][] temp;
		for (int i = 0; i < numSources; i++) {
			temp = ateCalcs[i].computeLocalFromPreviousObservations(states,
					cleanedOffsets[i][ROW_INDEX],
					cleanedOffsets[i][COLUMN_INDEX]);
			if (computeMultiInfoCoherence) {
				// Keep a link to this apparent transfer
				localComponents[1 + numSources] = temp;
			}
			try {
				local = MatrixUtils.add(local, temp);
			} catch (Exception e) {
				// Exception only thrown where arrays were not
				//  of the same length - should not happen here.
				return null;
			}
		}
		// Set statistics and if not periodic boundary conditions, 
		//  clean up points which don't get all information contributors
		setStatistics(local, cleanedOffsets, localComponents);
		localStatsValid = true;
		return local;
	}
	
	@Override
	public double[][] computeLocal(int states[][], int[] offsetOfDestFromSources) {
		double[][] local = aiCalc.computeLocal(states);
		double[][][] localComponents = null;
		if (computeMultiInfoCoherence) {
			localComponents = new double[numSources + 1][][];
			// Keep a link to the local active info
			localComponents[0] = local;
		}
		int[] cleanedOffsets = cleanOffsetOfDestFromSources(offsetOfDestFromSources);
		TransferEntropyCalculatorDiscrete ateCalc = 
			TransferEntropyCalculatorDiscrete.newInstance(base, k);
		ateCalc.setPeriodicBoundaryConditions(periodicBoundaryConditions);
		double[][] temp;
		for (int i = 0; i < numSources; i++) {
			temp = ateCalc.computeLocal(
					states, cleanedOffsets[i]);
			if (computeMultiInfoCoherence) {
				// Keep a link to this apparent transfer
				localComponents[1 + numSources] = temp;
			}
			try {
				local = MatrixUtils.add(local, temp);
			} catch (Exception e) {
				// Exception only thrown where arrays were not
				//  of the same length - should not happen here.
				return null;
			}
		}
		// Set statistics and if not periodic boundary conditions, 
		//  clean up points which don't get all information contributors
		setStatistics(local, cleanedOffsets, localComponents);
		localStatsValid = true;
		return local;
	}
	
	@Override
	public double[][][] computeLocal(int states[][][], int[][] offsetOfDestFromSources) {
		double[][][] local = aiCalc.computeLocal(states);
		double[][][][] localComponents = null;
		if (computeMultiInfoCoherence) {
			localComponents = new double[numSources + 1][][][];
			// Keep a link to the local active info
			localComponents[0] = local;
		}
		int[][] cleanedOffsets = cleanOffsetOfDestFromSources(offsetOfDestFromSources);
		TransferEntropyCalculatorDiscrete ateCalc = 
			TransferEntropyCalculatorDiscrete.newInstance(base, k);
		ateCalc.setPeriodicBoundaryConditions(periodicBoundaryConditions);
		double[][][] temp;
		for (int i = 0; i < numSources; i++) {
			temp = ateCalc.computeLocal(states,
					cleanedOffsets[i][ROW_INDEX],
					cleanedOffsets[i][COLUMN_INDEX]);
			if (computeMultiInfoCoherence) {
				// Keep a link to this apparent transfer
				localComponents[1 + numSources] = temp;
			}
			try {
				local = MatrixUtils.add(local, temp);
			} catch (Exception e) {
				// Exception only thrown where arrays were not
				//  of the same length - should not happen here.
				return null;
			}
		}
		// Set statistics and if not periodic boundary conditions, 
		//  clean up points which don't get all information contributors
		setStatistics(local, cleanedOffsets, localComponents);
		localStatsValid = true;
		return local;
	}
	
	@Override
	public double computeAverageLocal(int[][] states, int[] sourceOffsets) {
		// One could allow the call to defer here, however it will
		//  be much less memory intensive if we can use a single 
		//  ApparentTransferEntropy object
		average = aiCalc.computeAverageLocal(states);
		int[] cleanedOffsets = cleanOffsetOfDestFromSources(sourceOffsets);
		TransferEntropyCalculatorDiscrete ateCalc = 
			TransferEntropyCalculatorDiscrete.newInstance(base, k);
		ateCalc.setPeriodicBoundaryConditions(periodicBoundaryConditions);
		for (int i = 0; i < numSources; i++) {
			average += ateCalc.computeAverageLocal(states, cleanedOffsets[i]);
		}
		localStatsValid = false;
		return average;
	}
	
	@Override
	public double computeAverageLocal(int[][][] states, int[][] sourceOffsets) {
		// One could allow the call to defer here, however it will
		//  be much less memory intensive if we can use a single 
		//  ApparentTransferEntropy object
		average = aiCalc.computeAverageLocal(states);
		int[][] cleanedOffsets = cleanOffsetOfDestFromSources(sourceOffsets);
		TransferEntropyCalculatorDiscrete ateCalc = 
			TransferEntropyCalculatorDiscrete.newInstance(base, k);
		ateCalc.setPeriodicBoundaryConditions(periodicBoundaryConditions);
		for (int i = 0; i < numSources; i++) {
			average += ateCalc.computeAverageLocal(states,
					cleanedOffsets[i][ROW_INDEX],
					cleanedOffsets[i][COLUMN_INDEX]);
		}
		localStatsValid = false;
		return average;
	}
	
	@Override
	public double[] computeLocal(int states[][], int destCol, int[] sourcesAbsolute) {
		double[] local = aiCalc.computeLocal(states, destCol);
		double[][] localComponents = null;
		if (computeMultiInfoCoherence) {
			localComponents = new double[numSources + 1][];
			// Keep a link to the local active info
			localComponents[0] = local;
		}
		int[] cleanedOffsets = cleanAbsoluteSources(sourcesAbsolute, destCol);
		TransferEntropyCalculatorDiscrete ateCalc = 
			TransferEntropyCalculatorDiscrete.newInstance(base, k);
		ateCalc.setPeriodicBoundaryConditions(periodicBoundaryConditions);
		double[] temp;
		for (int i = 0; i < numSources; i++) {
			temp = ateCalc.computeLocal(
					states, cleanedOffsets[i], destCol);
			if (computeMultiInfoCoherence) {
				// Keep a link to this apparent transfer
				localComponents[1 + numSources] = temp;
			}
			try {
				local = MatrixUtils.add(local, temp);
			} catch (Exception e) {
				// Exception only thrown where arrays were not
				//  of the same length - should not happen here.
				return null;
			}
		}
		// Set statistics
		setStatistics(local, localComponents);
		localStatsValid = true;
		return local;
	}
	
	@Override
	public double[] computeLocal(int states[][][], int destAgentRow,
			int destAgentColumn, int[][] sourcesAbsolute) {
		double[] local = aiCalc.computeLocal(states, destAgentRow, destAgentColumn);
		double[][] localComponents = null;
		if (computeMultiInfoCoherence) {
			localComponents = new double[numSources + 1][];
			// Keep a link to the local active info
			localComponents[0] = local;
		}
		int[][] cleanedSourcesAbsolute = cleanAbsoluteSources(sourcesAbsolute, destAgentRow, destAgentColumn);
		TransferEntropyCalculatorDiscrete ateCalc = 
			TransferEntropyCalculatorDiscrete.newInstance(base, k);
		ateCalc.setPeriodicBoundaryConditions(periodicBoundaryConditions);
		double[] temp;
		for (int i = 0; i < numSources; i++) {
			temp = ateCalc.computeLocal(states, cleanedSourcesAbsolute[i][ROW_INDEX], cleanedSourcesAbsolute[i][COLUMN_INDEX],
					destAgentRow,
					destAgentColumn);
			if (computeMultiInfoCoherence) {
				// Keep a link to this apparent transfer
				localComponents[1 + numSources] = temp;
			}
			try {
				local = MatrixUtils.add(local, temp);
			} catch (Exception e) {
				// Exception only thrown where arrays were not
				//  of the same length - should not happen here.
				return null;
			}
		}
		// Set statistics
		setStatistics(local, localComponents);
		localStatsValid = true;
		return local;
	}

	@Override
	public double computeAverageLocal(int[][] states, int destCol, int[] sourcesAbsolute) {
		// One could allow the call to defer here, however it will
		//  be much less memory intensive if we can use a single 
		//  ApparentTransferEntropy object
		average = aiCalc.computeAverageLocal(states, destCol);
		int[] cleanedSourcesAbsolute = cleanAbsoluteSources(sourcesAbsolute, destCol);
		TransferEntropyCalculatorDiscrete ateCalc = 
			TransferEntropyCalculatorDiscrete.newInstance(base, k);
		ateCalc.setPeriodicBoundaryConditions(periodicBoundaryConditions);
		for (int i = 0; i < numSources; i++) {
			average += ateCalc.computeAverageLocal(states, cleanedSourcesAbsolute[i], destCol);
		}
		localStatsValid = false;
		return average;
	}

	@Override
	public double computeAverageLocal(int states[][][], int destAgentRow,
			int destAgentColumn, int[][] sourcesAbsolute) {
		// One could allow the call to defer here, however it will
		//  be much less memory intensive if we can use a single 
		//  ApparentTransferEntropy object
		average = aiCalc.computeAverageLocal(states);
		int[][] cleanedSourcesAbsolute = cleanAbsoluteSources(sourcesAbsolute, destAgentRow, destAgentColumn);
		TransferEntropyCalculatorDiscrete ateCalc = 
			TransferEntropyCalculatorDiscrete.newInstance(base, k);
		ateCalc.setPeriodicBoundaryConditions(periodicBoundaryConditions);
		for (int i = 0; i < numSources; i++) {
			average += ateCalc.computeAverageLocal(states,
					cleanedSourcesAbsolute[i][ROW_INDEX],
					cleanedSourcesAbsolute[i][COLUMN_INDEX]);
		}
		localStatsValid = false;
		return average;
	}
	
	@Override
	public void setPeriodicBoundaryConditions(boolean periodicBoundaryConditions) {
		super.setPeriodicBoundaryConditions(periodicBoundaryConditions);
		if (ateCalcs != null) {
			for (int i = 0; i < numSources; i++) {
				ateCalcs[i].setPeriodicBoundaryConditions(periodicBoundaryConditions);
			}
		}
	}

	@Override
	public double getLastMax() {
		if (localStatsValid) {
			return super.getLastMax();
		}
		throw new RuntimeException("Last maximum is not valid after previous computation by SeparableInfoCalculatorByAddition");
	}

	@Override
	public double getLastMin() {
		if (localStatsValid) {
			return super.getLastMin();
		}
		throw new RuntimeException("Last minimum is not valid after previous computation by SeparableInfoCalculatorByAddition");
	}

	@Override
	public double getLastAverageNegative() {
		if (localStatsValid) {
			return super.getLastAverageNegative();
		}
		throw new RuntimeException("Last average negative sep is not valid after previous computation by SeparableInfoCalculatorByAddition");
	}

	@Override
	public double getLastAveragePositive() {
		if (localStatsValid) {
			return super.getLastAveragePositive();
		}
		throw new RuntimeException("Last average positive sep is not valid after previous computation by SeparableInfoCalculatorByAddition");
	}

	private void setStatistics(double[] localValues, double[][] localComponents) {
		int timeSteps = localValues.length;

		// Compute average, max, min, avPositiveLocal, avNegativeLocal
		average = 0;
		max = localValues[k];
		min = localValues[k];
		avPositiveLocal = 0;
		avNegativeLocal = 0;
		for (int t = k; t < timeSteps; t++) {
			average += localValues[t];
			if (localValues[t] > 0) {
				avPositiveLocal += localValues[t];
			} else {
				avNegativeLocal += localValues[t];
			}
			if (localValues[t] > max) {
				max = localValues[t];
			} else if (localValues[t] < min) {
				min = localValues[t];
			}
		}
		average /= (double) (timeSteps - k);
		avPositiveLocal /= (double) (timeSteps - k);
		avNegativeLocal /= (double) (timeSteps - k);
		
		// Now compute the coherence of computation
		if (computeMultiInfoCoherence) {
			miCalc.startAddObservations();
			double[] miTuple = new double[numSources + 1];
			for (int t = k; t < timeSteps; t++) {
				// Construct the multi-info tuple
				for (int i = 0; i < numSources + 1; i++) {
					miTuple[i] = localComponents[i][t];
				}
				// Add this tuple to the observations
				miCalc.addObservation(miTuple);
			}
			try {
				miCalc.finaliseAddObservations();
			} catch (Exception e) {
				// an exception would only be thrown if we changed the number of causal contributors here
				//  which simply will not happen. Just in case it does, we'll throw a runtime exception
				throw new RuntimeException("Number of causal contributors changed from intialisation to calculation!");
			}
		}
	}

	private void setStatistics(double[][] localValues, int[] cleanedSourcesOffsets,
			double[][][] localComponents) {
		int timeSteps = localValues.length;
		int numAgents = localValues[0].length;
		
		// Clean up the point which didn't have all of the local information contributors
		// if we're doing non-peridoic boundary conditions
		int minAgentOffset = MatrixUtils.min(cleanedSourcesOffsets);
		int maxAgentOffset = MatrixUtils.max(cleanedSourcesOffsets);
		int nonPeriodicStartAgent = Math.max(0, maxAgentOffset);
		int nonPeriodicEndAgent = numAgents - 1 + Math.min(0, minAgentOffset);

		// Clean up if required
		if (!periodicBoundaryConditions) {
			for (int t = k; t < timeSteps; t++) {
				for (int r = 0; r < nonPeriodicStartAgent; r++) {
					localValues[t][r] = 0;
				}
				for (int r = nonPeriodicEndAgent + 1; r < numAgents; r++) {
					localValues[t][r] = 0;
				}
			}
		}

		// Compute average, max, min, avPositiveLocal, avNegativeLocal
		average = 0;
		max = localValues[k][0];
		min = localValues[k][0];
		avPositiveLocal = 0;
		avNegativeLocal = 0;
		for (int t = k; t < timeSteps; t++) {
			for (int r = periodicBoundaryConditions ? 0 : nonPeriodicStartAgent;
			 r < (periodicBoundaryConditions ? numAgents : nonPeriodicEndAgent + 1);
			 r++) {
				average += localValues[t][r];
				if (localValues[t][r] > 0) {
					avPositiveLocal += localValues[t][r];
				} else {
					avNegativeLocal += localValues[t][r];
				}
				if (localValues[t][r] > max) {
					max = localValues[t][r];
				} else if (localValues[t][r] < min) {
					min = localValues[t][r];
				}
			}
		}
		if (periodicBoundaryConditions) {
			average /= (double) ((timeSteps - k) * numAgents);
			avPositiveLocal /= (double) ((timeSteps - k) * numAgents);
			avNegativeLocal /= (double) ((timeSteps - k) * numAgents);
		} else {
			average /= (double) ((timeSteps - k) * (nonPeriodicEndAgent - nonPeriodicStartAgent + 1));
			avPositiveLocal /= (double) ((timeSteps - k) * (nonPeriodicEndAgent - nonPeriodicStartAgent + 1));
			avNegativeLocal /= (double) ((timeSteps - k) * (nonPeriodicEndAgent - nonPeriodicStartAgent + 1));
		}

		// Now compute the coherence of computation
		if (computeMultiInfoCoherence) {
			miCalc.startAddObservations();
			double[] miTuple = new double[numSources + 1];
			for (int t = k; t < timeSteps; t++) {
				for (int r = periodicBoundaryConditions ? 0 : nonPeriodicStartAgent;
				 r < (periodicBoundaryConditions ? numAgents : nonPeriodicEndAgent + 1);
				 r++) {
					// Construct the multi-info tuple
					for (int i = 0; i < numSources + 1; i++) {
						miTuple[i] = localComponents[i][t][r];
					}
					// Add this tuple to the observations
					miCalc.addObservation(miTuple);
				}
			}
			try {
				miCalc.finaliseAddObservations();
			} catch (Exception e) {
				// an exception would only be thrown if we changed the number of causal contributors here
				//  which simply will not happen. Just in case it does, we'll throw a runtime exception
				throw new RuntimeException("Number of causal contributors changed from intialisation to calculation!");
			}
		}
	}

	private void setStatistics(double[][][] localValues, int[][] cleanedSourcesOffsets,
			double[][][][] localComponents) {
		int timeSteps = localValues.length;
		int numAgentRows = localValues[0].length;
		int numAgentColumns = localValues[0][0].length;
		
		// Clean up the point which didn't have all of the local information contributors
		// if we're doing non-peridoic boundary conditions
		int minRowOffset = MatrixUtils.min(cleanedSourcesOffsets, ROW_INDEX);
		int maxRowOffset = MatrixUtils.max(cleanedSourcesOffsets, ROW_INDEX);
		int nonPeriodicStartRow = Math.max(0, maxRowOffset);
		int nonPeriodicEndRow = numAgentRows - 1 + Math.min(0, minRowOffset);
		int minColumnOffset = MatrixUtils.min(cleanedSourcesOffsets, COLUMN_INDEX);
		int maxColumnOffset = MatrixUtils.max(cleanedSourcesOffsets, COLUMN_INDEX);
		int nonPeriodicStartColumn = Math.max(0, maxColumnOffset);
		int nonPeriodicEndColumn = numAgentColumns - 1 + Math.min(0, minColumnOffset);

		System.out.println(periodicBoundaryConditions + " " + nonPeriodicStartRow + " " +
				nonPeriodicEndRow + " " + nonPeriodicStartColumn + " " + nonPeriodicEndColumn);
		
		// Clean up if required
		if (!periodicBoundaryConditions) {
			for (int t = k; t < timeSteps; t++) {
				for (int r = 0; r < nonPeriodicStartRow; r++) {
					for (int c = 0; c < nonPeriodicStartColumn; c++) {
						localValues[t][r][c] = 0;
					}
					for (int c = nonPeriodicEndColumn + 1; c < numAgentColumns; c++) {
						localValues[t][r][c] = 0;
					}
				}
				for (int r = nonPeriodicEndRow + 1; r < numAgentRows; r++) {
					for (int c = 0; c < nonPeriodicStartColumn; c++) {
						localValues[t][r][c] = 0;
					}
					for (int c = nonPeriodicEndColumn + 1; c < numAgentColumns; c++) {
						localValues[t][r][c] = 0;
					}
				}
			}
		}
		
		// Compute average, max, min, avPositiveLocal, avNegativeLocal
		average = 0;
		max = localValues[k][0][0];
		min = localValues[k][0][0];
		avPositiveLocal = 0;
		avNegativeLocal = 0;
		for (int t = k; t < timeSteps; t++) {
			for (int r = periodicBoundaryConditions ? 0 : nonPeriodicStartRow;
			 r < (periodicBoundaryConditions ? numAgentRows : nonPeriodicEndRow + 1);
			 r++) {
				for (int c = periodicBoundaryConditions ? 0 : nonPeriodicStartColumn;
				c < (periodicBoundaryConditions ? numAgentColumns : nonPeriodicEndColumn + 1);
				c++) {
					average += localValues[t][r][c];
					if (localValues[t][r][c] > 0) {
						avPositiveLocal += localValues[t][r][c];
					} else {
						avNegativeLocal += localValues[t][r][c];
					}
					if (localValues[t][r][c] > max) {
						max = localValues[t][r][c];
					} else if (localValues[t][r][c] < min) {
						min = localValues[t][r][c];
					}
				}
			}
		}
		if (periodicBoundaryConditions) {
			average /= (double) ((timeSteps - k) * numAgentRows * numAgentColumns);
			avPositiveLocal /= (double) ((timeSteps - k) * numAgentRows * numAgentColumns);
			avNegativeLocal /= (double) ((timeSteps - k) * numAgentRows * numAgentColumns);
		} else {
			average /= (double) ((timeSteps - k) * (nonPeriodicEndRow - nonPeriodicStartRow + 1) *
					(nonPeriodicEndColumn - nonPeriodicStartColumn + 1));
			avPositiveLocal /= (double) ((timeSteps - k) * (nonPeriodicEndRow - nonPeriodicStartRow + 1) *
					(nonPeriodicEndColumn - nonPeriodicStartColumn + 1));
			avNegativeLocal /= (double) ((timeSteps - k) * (nonPeriodicEndRow - nonPeriodicStartRow + 1) *
					(nonPeriodicEndColumn - nonPeriodicStartColumn + 1));
		}
		
		// Now compute the coherence of computation
		if (computeMultiInfoCoherence) {
			miCalc.startAddObservations();
			double[] miTuple = new double[numSources + 1];
			for (int t = k; t < timeSteps; t++) {
				for (int r = periodicBoundaryConditions ? 0 : nonPeriodicStartRow;
				 r < (periodicBoundaryConditions ? numAgentRows : nonPeriodicEndRow + 1);
				 r++) {
					for (int c = periodicBoundaryConditions ? 0 : nonPeriodicStartColumn;
					c < (periodicBoundaryConditions ? numAgentColumns : nonPeriodicEndColumn + 1);
					c++) {
						// Construct the multi-info tuple
						for (int i = 0; i < numSources + 1; i++) {
							miTuple[i] = localComponents[i][t][r][c];
						}
						// Add this tuple to the observations
						miCalc.addObservation(miTuple);
					}
				}
			}
			try {
				miCalc.finaliseAddObservations();
			} catch (Exception e) {
				// an exception would only be thrown if we changed the number of causal contributors here
				//  which simply will not happen. Just in case it does, we'll throw a runtime exception
				throw new RuntimeException("Number of causal contributors changed from intialisation to calculation!");
			}
		}
	}

	@Override
	public boolean canComputeMultiInfoCoherenceFromAverageOfObservations() {
		return false;
	}
	
	
}
