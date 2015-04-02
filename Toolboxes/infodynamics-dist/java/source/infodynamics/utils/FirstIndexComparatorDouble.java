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
 * <p>Class used to sort an array of array of doubles based on the first element in each array.
 * </p> 
 * 
 * <p>Usage: <code>java.util.Arrays.sort(double[][], FirstIndexComparatorDouble.getInstance())</code>;
 * </p>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class FirstIndexComparatorDouble implements java.util.Comparator<double[]> {
	
	private static FirstIndexComparatorDouble instance = null;
	
	private FirstIndexComparatorDouble() {
	}

	public static FirstIndexComparatorDouble getInstance() {
		if (instance == null) {
			instance = new FirstIndexComparatorDouble();
		}
		return instance;
	}
	
	public int compare(double[] a1, double[] a2) {
		if (a1[0] < a2[0]) {
			return -1;
		} else if (a1[0] == a2[0]) {
			return 0;
		} else {
			return 1;
		}
	}
}