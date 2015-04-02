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

/**
 * A base class for calculators computing measures for
 * a single variable which
 * require knowledge of the embedded past state of a univariate
 * discrete (ie int[]) variable.
 * 
 * <p>This combines functionality for single agents from
 * {@link SingleAgentMeasureDiscrete} with functionality
 * required in the context of the past provided by
 * {@link ContextOfPastMeasureCalculatorDiscrete}.</p>
 * 
 * <p>Usage is as defined in {@link InfoMeasureCalculatorDiscrete}, with
 * extra methods for supplying observations and making 
 * calculations defined in {@link SingleAgentMeasureDiscrete}</p>.
 * 
 * <p>Users should not need to deal with this class directly;
 * it is simply used to gather common functionality for several
 * child classes.
 * </p>
 * 
 * TODO Make the Active info storage and entropy calculators inherit from this
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public abstract class SingleAgentMeasureDiscreteInContextOfPastCalculator extends
		ContextOfPastMeasureCalculatorDiscrete implements SingleAgentMeasureDiscrete {

	/**
	 * Construct the calculator
	 * 
	 * @param base number of quantisation levels for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param history embedding length
	 */
	public SingleAgentMeasureDiscreteInContextOfPastCalculator(int base, int history) {
		super(base, history);
	}
	
	@Override
	public final double[] computeLocal(int[] states) {
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}

	@Override
	public final double[][] computeLocal(int[][] states) {
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}

	@Override
	public final double[][][] computeLocal(int[][][] states) {
		initialise();
		addObservations(states);
		return computeLocalFromPreviousObservations(states);
	}

	@Override
	public final double computeAverageLocal(int[] states) {
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}	

	@Override
	public final double computeAverageLocal(int[][] states) {
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}

	@Override
	public final double computeAverageLocal(int[][][] states) {
		initialise();
		addObservations(states);
		return computeAverageLocalOfObservations();
	}

	@Override
	public final double[] computeLocal(int[][] states, int col) {
		initialise();
		addObservations(states, col);
		return computeLocalFromPreviousObservations(states, col);
	}

	@Override
	public final double[] computeLocal(int[][][] states, int index1, int index2) {
		initialise();
		addObservations(states, index1, index2);
		return computeLocalFromPreviousObservations(states, index1, index2);
	}

	@Override
	public final double computeAverageLocal(int[][] states, int col) {
		initialise();
		addObservations(states, col);
		return computeAverageLocalOfObservations();
	}
	
	@Override
	public final double computeAverageLocal(int[][][] states, int index1, int index2) {
		initialise();
		addObservations(states, index1, index2);
		return computeAverageLocalOfObservations();
	}
}
