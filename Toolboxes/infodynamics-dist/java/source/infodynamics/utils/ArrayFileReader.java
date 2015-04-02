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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/**
 * Reads basic files with one or two dimensional arrays of integers
 *  and doubles,
 *  delimited by space, tab or comma.
 * One dimensional arrays should have a single number on each line.
 * Files can have comment lines, which start with # or % characters.
 * Empty lines are allowed and skipped.
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class ArrayFileReader {
	private String filename;
	private BufferedReader br;
	private int rows;
	private int columns;
	private boolean isOpen = false;
	
	static final String DELIMITER = "[ \t,]";
	
	/**
	 * Creates the ArrayFileReader and opens the given file for reading.
	 * 
	 * @param arrayFilename
	 */
	public ArrayFileReader(String arrayFilename) {
		filename = arrayFilename;
	}
	
	/**
	 * Creates the ArrayFileReader with a {@link File} instance
	 * 
	 * @param file {@link File} instance to read
	 */
	public ArrayFileReader(File file) {
		filename = file.getPath();
	}
	
	private void openArrayFile()
		throws FileNotFoundException, IOException {
		if (isOpen) {
			return;
		}
		
		// Work out how many columns and how many rows we have
		br = new BufferedReader(new FileReader(filename));
		rows = 0;
		columns = 0;
		for (String line = br.readLine(); line != null; line = br.readLine()) {
			// Strip leading whitespace
			line = line.replaceFirst("^[ \t]*", "");
			if (line.startsWith("#") || line.startsWith("%")) {
				// comment line
				continue;
			}
			String[] stringValues = line.split(DELIMITER);
			if (stringValues.length == 0) {
				// Nothing on this line
				continue;
			}
			// Check that there is something other than whitespace here
			int columnsHere = 0;
			for (int i = 0; i < stringValues.length; i++) {
				if (!stringValues[i].equals("")) {
					// Got a non-empty entity
					columnsHere++;
				}
			}
			if (columnsHere == 0) {
				// Nothing on this line
				continue;
			}
			// Got something here
			if (columns == 0) {
				// Count the number of columns on the first line
				columns = columnsHere;
			}
			// Now make sure that the number of columns here
			//  matches the number on previous lines
			if (columnsHere != columns) {
				throw new IOException("Number of columns " + columnsHere +
						" on row " + rows +
						" does not match number on previous lines " + columns);
			}
			// Else this row is OK
			rows++;
		}
		br.close();
		
		// Now have the file open and ready to read in
		br = new BufferedReader(new FileReader(filename));
		isOpen = true;
	}
	
	private void closeArrayFile() throws IOException {
		br.close();
		isOpen = false;
	}

	/**
	 * Attempt to read a 2D int matrix from the file
	 * 
	 * @return the 2D int[][] matrix
	 * @throws Exception
	 */
	public int[][] getInt2DMatrix() throws Exception {
		openArrayFile();
		// Create the return array
		int[][] values = new int[rows][columns];
		
		int r = 0;
		for (String line = br.readLine(); line != null; line = br.readLine()) {
			// Strip leading whitespace
			line = line.replaceFirst("^[ \t]*", "");
			if (line.startsWith("#") || line.startsWith("%")) {
				// comment line
				continue;
			}
			String[] stringValues = line.split(DELIMITER);
			if (stringValues.length == 0) {
				// Nothing on this line
				continue;
			}
			// Check that there is something other than whitespace here
			int c = 0;
			for (int i = 0; i < stringValues.length; i++) {
				if (!stringValues[i].equals("")) {
					// Got a non-empty entity
					// Not checking array bounds since this
					//  was already checked by openArrayFile().
					values[r][c] = Integer.parseInt(stringValues[i]);
					c++;
				}
			}
			if (c == 0) {
				// Nothing was on this line
				continue;
			}
			// Got some values here.
			// Now make sure that the number of columns here
			//  matches the number on previous lines
			if (c != columns) {
				// This should not occur as it was already checked
				//  in the openArrayFile() method
				throw new IOException("Number of columns " + c +
						" on row " + rows +
						" does not match number on previous lines " + columns);
			}
			// Else this row is OK
			r++;
		}

		closeArrayFile();
		return values;
	}

	/**
	 * Attempt to read a 1D int array from the file
	 * 
	 * @return the array that was read
	 * @throws Exception
	 */
	public int[] getIntArray() throws Exception {
		openArrayFile();
		
		if ((rows != 1) && (columns != 1)) {
			throw new Exception("Cannot get an int array when neither the number of rows nor columns equals 1");
		}
		boolean alongRows = (rows > columns);
		int length = alongRows ? rows : columns;
		
		// Create the return array
		int[] values = new int[length];
		
		if (alongRows) {
			int r = 0;
			for (String line = br.readLine(); line != null; line = br.readLine()) {
				// Strip leading whitespace
				line = line.replaceFirst("^[ \t]*", "");
				if (line.startsWith("#") || line.startsWith("%")) {
					// comment line
					continue;
				}
				String[] stringValues = line.split(DELIMITER);
				if (stringValues.length == 0) {
					// Nothing on this line
					continue;
				}
				// Check that there is something other than whitespace here
				boolean readOneOnThisLine = false;
				for (int i = 0; i < stringValues.length; i++) {
					if (!stringValues[i].equals("")) {
						if (readOneOnThisLine) {
							// We've already read an integer on this line - this is an error condition
							throw new IOException("Cannot have multiple values on each line for a column vector");
						}
						// Got a non-empty entity
						// Not checking array bounds since this
						//  was already checked by openArrayFile().
						values[r] = Integer.parseInt(stringValues[i]);
						readOneOnThisLine = true;
					}
				}
				if (!readOneOnThisLine) {
					// Nothing was on this line
					continue;
				}
				// Got some value here.
				r++;
			}
		} else {
			boolean readFromOneLine = false;
			for (String line = br.readLine(); line != null; line = br.readLine()) {
				// Strip leading whitespace
				line = line.replaceFirst("^[ \t]*", "");
				if (line.startsWith("#") || line.startsWith("%")) {
					// comment line
					continue;
				}
				String[] stringValues = line.split(DELIMITER);
				if (stringValues.length == 0) {
					// Nothing on this line
					continue;
				}
				// Check that there is something other than whitespace here
				int c = 0;
				for (int i = 0; i < stringValues.length; i++) {
					if (!stringValues[i].equals("")) {
						if ((c == 0) && readFromOneLine) {
							// We're about to start reading from this line, but it's not the first
							//  line we've read from
							throw new IOException("Cannot have data on multiple lines for a row vector");
						}
						// Got a non-empty entity - read all from this line
						// Not checking array bounds since this
						//  was already checked by openArrayFile().
						values[c] = Integer.parseInt(stringValues[i]);
						readFromOneLine = true;
						c++;
					}
				}
				if (c == 0) {
					// Nothing was on this line
					continue;
				}
				// Got at least one value here - check if we read enough
				if (c != columns) {
					// We did not read the correct number of columns here:
					throw new IOException("Not enough elements in the row vector");
				}
			}
		}
		
		closeArrayFile();
		return values;
	}

	/**
	 * Attempt to read a double[][] array from the file
	 * 
	 * @return the read double[][] array
	 * @throws Exception
	 */
	public double[][] getDouble2DMatrix() throws Exception {
		openArrayFile();
		// Create the return array
		double[][] values = new double[rows][columns];
		
		int r = 0;
		for (String line = br.readLine(); line != null; line = br.readLine()) {
			// Strip leading whitespace
			line = line.replaceFirst("^[ \t]*", "");
			if (line.startsWith("#") || line.startsWith("%")) {
				// comment line
				continue;
			}
			String[] stringValues = line.split(DELIMITER);
			if (stringValues.length == 0) {
				// Nothing on this line
				continue;
			}
			// Check that there is something other than whitespace here
			int c = 0;
			for (int i = 0; i < stringValues.length; i++) {
				if (!stringValues[i].equals("")) {
					// Got a non-empty entity
					// Not checking array bounds since this
					//  was already checked by openArrayFile().
					values[r][c] = Double.parseDouble(stringValues[i]);
					c++;
				}
			}
			if (c == 0) {
				// Nothing was on this line
				continue;
			}
			// Got some values here.
			// Now make sure that the number of columns here
			//  matches the number on previous lines
			if (c != columns) {
				// This should not occur as it was already checked
				//  in the openArrayFile() method
				throw new IOException("Number of columns " + c +
						" on row " + rows +
						" does not match number on previous lines " + columns);
			}
			// Else this row is OK
			r++;
		}

		closeArrayFile();
		return values;
	}

	/**
	 * Attempt to read a double[] array from the file
	 * 
	 * @return the read double[] array
	 * @throws Exception
	 */
	public double[] getDoubleArray() throws Exception {
		openArrayFile();
		
		if ((rows != 1) && (columns != 1)) {
			throw new Exception("Cannot get an double array when neither the number of rows nor columns equals 1");
		}
		boolean alongRows = (rows > columns);
		int length = alongRows ? rows : columns;
		
		// Create the return array
		double[] values = new double[length];
		
		if (alongRows) {
			int r = 0;
			for (String line = br.readLine(); line != null; line = br.readLine()) {
				// Strip leading whitespace
				line = line.replaceFirst("^[ \t]*", "");
				if (line.startsWith("#") || line.startsWith("%")) {
					// comment line
					continue;
				}
				String[] stringValues = line.split(DELIMITER);
				if (stringValues.length == 0) {
					// Nothing on this line
					continue;
				}
				// Check that there is something other than whitespace here
				boolean readOneOnThisLine = false;
				for (int i = 0; i < stringValues.length; i++) {
					if (!stringValues[i].equals("")) {
						if (readOneOnThisLine) {
							// We've already read an double on this line - this is an error condition
							throw new IOException("Cannot have multiple values on each line for a column vector");
						}
						// Got a non-empty entity
						// Not checking array bounds since this
						//  was already checked by openArrayFile().
						values[r] = Double.parseDouble(stringValues[i]);
						readOneOnThisLine = true;
					}
				}
				if (!readOneOnThisLine) {
					// Nothing was on this line
					continue;
				}
				// Got some value here.
				r++;
			}
		} else {
			boolean readFromOneLine = false;
			for (String line = br.readLine(); line != null; line = br.readLine()) {
				// Strip leading whitespace
				line = line.replaceFirst("^[ \t]*", "");
				if (line.startsWith("#") || line.startsWith("%")) {
					// comment line
					continue;
				}
				String[] stringValues = line.split(DELIMITER);
				if (stringValues.length == 0) {
					// Nothing on this line
					continue;
				}
				// Check that there is something other than whitespace here
				int c = 0;
				for (int i = 0; i < stringValues.length; i++) {
					if (!stringValues[i].equals("")) {
						if ((c == 0) && readFromOneLine) {
							// We're about to start reading from this line, but it's not the first
							//  line we've read from
							throw new IOException("Cannot have data on multiple lines for a row vector");
						}
						// Got a non-empty entity - read all from this line
						// Not checking array bounds since this
						//  was already checked by openArrayFile().
						values[c] = Double.parseDouble(stringValues[i]);
						readFromOneLine = true;
						c++;
					}
				}
				if (c == 0) {
					// Nothing was on this line
					continue;
				}
				// Got at least one value here - check if we read enough
				if (c != columns) {
					// We did not read the correct number of columns here:
					throw new IOException("Not enough elements in the row vector");
				}
			}
		}
		
		closeArrayFile();
		return values;
	}

}
