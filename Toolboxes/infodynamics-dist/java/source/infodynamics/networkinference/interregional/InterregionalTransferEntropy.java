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

import infodynamics.measures.continuous.TransferEntropyCalculatorMultiVariate;
import infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorMultiVariateKraskov;
import infodynamics.utils.ParsedProperties;

import java.util.Vector;


/**
 * <p>Compute a significance score for the inter-regional transfer entropy
 * between two sets, based on looking at transfer between elements taken 
 * e at a time.
 * </p> 
 * 
 * @author Joseph Lizier, joseph.lizier at gmail.com 
 *
 */
public class InterregionalTransferEntropy extends InterregionalChannelMeasure {

	protected int k;
	protected static final String PROP_K_TE = "props.interregionalChannel.te.k";

	public InterregionalTransferEntropy() {
		super();
	}
	
	public void initialise(ParsedProperties props) throws Exception {
		setK(props.getIntProperty(PROP_K_TE));
		super.initialise(props);
	}

	public void initialise(int k, int jointVars1, int jointVars2, int maxNumSubsets) throws Exception {
		setK(k);
		super.initialise(jointVars1, jointVars2, maxNumSubsets);
	}

	protected void initialiseCalculator() throws Exception {
		TransferEntropyCalculatorMultiVariate teChannelCalc =
			(TransferEntropyCalculatorMultiVariate) channelCalc;
		teChannelCalc.initialise(k, jointVars1, jointVars2);
	}

	public void setK(int k) {
		this.k = k;
	}

	@Override
	public int[] computeTimeIndicesForLocalValues() throws Exception {
		int[] localtimeIndices;
		if (allValid) {
			localtimeIndices = new int[region1.length - k];
			for (int t = k; t < region1.length; t++) {
				localtimeIndices[t - k] = t;
			}
		} else if (!validityForIndividualElements) {
			int numObs = computeNumberOfObservations(jointValidity1, jointValidity2);
			localtimeIndices = new int[numObs];
			TransferEntropyCalculatorMultiVariateKraskov tecmvKras =
				new TransferEntropyCalculatorMultiVariateKraskov();
			tecmvKras.initialise(k);
			Vector<int[]> startAndEndTimePairs = tecmvKras.computeStartAndEndTimePairs(
					jointValidity1, jointValidity2);
			// We've found the set of start and end times for this pair
			int obsNum = 0;
			for (int[] timePair : startAndEndTimePairs) {
				int startTime = timePair[0];
				int endTime = timePair[1];
				for (int t = startTime + k; t <= endTime; t++) {
					localtimeIndices[obsNum++] = t;
				}
			}

		} else {
			// TODO Need to implement this properly in the super class.
			// for the moment I am assuming it gets implemented by putting locals
			//  at all time points
			localtimeIndices = new int[region1.length - k];
			for (int t = k; t < region1.length; t++) {
				localtimeIndices[t - k] = t;
			}
		}
		return localtimeIndices;
	}	
	
	@Override
	protected int computeNumObservationsToReorder() throws Exception {
		if (allValid) {
			// all the observations are valid
			return region1.length - k;
		} else if (!validityForIndividualElements) {
			// we've got joint validity for each time series
			return computeNumberOfObservations(jointValidity1, jointValidity2);
		} else {
			// We've been given validity for each individual sub-variable.
			// There is a different number of observations for each sub-variable.
			// It is best if we just reorder all of the observations here
			//  for every subset.
			return region1.length;
		}
	}

	/**
	 * Dummy method to compute the number of observations for the given validities
	 * 
	 * @param sourceValid
	 * @param destValid
	 * @return
	 * @throws Exception
	 */
	private int computeNumberOfObservations(boolean[] sourceValid, boolean[] destValid) throws Exception {
		TransferEntropyCalculatorMultiVariateKraskov tecmvKras =
			new TransferEntropyCalculatorMultiVariateKraskov();
		
		tecmvKras.initialise(k);
		Vector<int[]> startAndEndTimePairs = tecmvKras.computeStartAndEndTimePairs(
				sourceValid, destValid);
		
		// We've found the set of start and end times for this pair
		int numObservations = 0;
		for (int[] timePair : startAndEndTimePairs) {
			int startTime = timePair[0];
			int endTime = timePair[1];
			numObservations += endTime - startTime + 1 - k;
		}
		return numObservations;
	}
}
