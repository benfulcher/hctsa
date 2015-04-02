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

import java.util.Arrays;

/**
 * IntArrayWrapper is used to wrap a long array so as to be able to give it a hash value
 *  for placement into a hash table.
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class IntArrayWrapper {

	private int[] array;
	private int firstCols = 0;
	
	public IntArrayWrapper() {
		super();
	}

	public IntArrayWrapper(int[] newArray) {
		super();
		setArray(newArray, newArray.length);
	}

	public IntArrayWrapper(int[] newArray, int firstcolumns) {
		super();
		setArray(newArray, firstcolumns);
	}

	public int hashCode() {
		return Arrays.hashCode(array);
	}
	
	public int[] getArray() {
		return array;
	}

	public void setArray(int[] newArray, int firstcolumns) {
		array = newArray;
		firstCols = firstcolumns;
	}
	
	public int getArrayUsedLength() {
		return firstCols;
	}
	
	public boolean equals(Object o) {
		if (o == this) {
			return true;
		}
		
		IntArrayWrapper iaw2 = (IntArrayWrapper) o;
		if (iaw2.getArrayUsedLength() != firstCols) {
			return false;
		}
		
		// check deeply inside arrays
		int[] inArray = iaw2.getArray();
		for (int i = 0; i < firstCols; i++) {
			if (array[i] != inArray[i]) {
				return false;
			}
		}
		return true;
	}
}
