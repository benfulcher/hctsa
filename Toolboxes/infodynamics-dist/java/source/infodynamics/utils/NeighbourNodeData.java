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

package infodynamics.utils;

/**
 * Class for storing nearest neighbour values during the search
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class NeighbourNodeData implements Comparable<NeighbourNodeData> {
	public int sampleIndex;
	public double[] norms; // norms in each high-level variable
	public double distance; // Assertion: distance is the max of norms
	
	/**
	 * Create an instance representing data about one given nearest neighbour
	 *  to another data point.
	 *  
	 * 
	 * @param sampleIndex index of the neighbour
	 * @param norms norms between the neighbour and the other data point,
	 *  for each high-level variable. 
	 * @param distance the max of the norms (used for sorting 
	 *  NeighbourNodeData objects in a PriorityQueue)
	 */
	public NeighbourNodeData(int sampleIndex, double[] norms, double distance) {
		super();
		this.norms = norms;
		this.sampleIndex = sampleIndex;
		this.distance = distance;
	}

	/**
	 * Override's {@link Comparable#compareTo(Object)} to provide
	 *  a natural comparison for the NeighbourNodeData class,
	 *  based on the underlying distance member.
	 *  
	 * <p><b>IMPORTANT</b> -- we reverse the usual comparison return
	 *  values, returning a positive number if other is greater,
	 *  and a negative number if other is less. This is so that a 
	 *  {@link java.util.PriorityQueue} of NeighbourNodeData objects
	 *  will hold that with the largest distance as the head (whereas
	 *  for standard return values from this method it would be
	 *  the other way around).
	 * 
	 * @param other Other NeighbourNodeData to compare to 
	 */
	@Override
	public int compareTo(NeighbourNodeData other) {
		if (distance < other.distance) {
			// Normally would return -1 here but we flip it -- see
			//  header comments
			return 1;
		} else if (distance > other.distance) {
			// Normally would return +1 here but we flip it -- see
			//  header comments
			return -1;
		}
		// distances are equal
		return 0;
	}
}