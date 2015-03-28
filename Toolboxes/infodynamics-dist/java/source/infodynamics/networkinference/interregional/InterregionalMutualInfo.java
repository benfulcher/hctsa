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

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.utils.ParsedProperties;

/**
 * <p>Compute a significance score for the inter-regional mutual information
 * between two sets, based on looking at transfer between elements taken 
 * e at a time.
 * </p> 
 * 
 * @author Joseph Lizier, joseph.lizier at gmail.com 
 *
 */
public class InterregionalMutualInfo extends InterregionalChannelMeasure {

	private int timeDiff = 0;
	
	public InterregionalMutualInfo() {
		super();
	}
	
	public void initialise(ParsedProperties props) throws Exception {
		String timeDiffPropName = InterregionalChannelMeasure.PROP_CALCULATOR_PROPERTIES_PREFIX + 
						MutualInfoCalculatorMultiVariate.PROP_TIME_DIFF;
		if (props.containsProperty(timeDiffPropName)) {
			timeDiff = props.getIntProperty(timeDiffPropName);
		}
		super.initialise(props);
	}

	@Override
	protected int computeNumObservationsToReorder() throws Exception {
		if (allValid) {
			// all the observations are valid
			return region1.length - timeDiff;
		} else if (!validityForIndividualElements) {
			// we've got joint validity for each time series
			// Check how many of the time steps are both valid
			int numValidTimePoints = 0;
			for (int t = 0; t < jointValidity1.length - timeDiff; t++) {
				if (jointValidity1[t] && jointValidity2[t + timeDiff]) {
					numValidTimePoints++;
				}
			}
			return numValidTimePoints;
		} else {
			// We've been given validity for each individual sub-variable.
			// There is a different number of observations for each sub-variable.
			// It is best if we just reorder all of the observations here
			//  for every subset.
			return region1.length;
		}
	}

	@Override
	public int[] computeTimeIndicesForLocalValues() throws Exception {
		throw new Exception("Not implemented yet");
	}
}
