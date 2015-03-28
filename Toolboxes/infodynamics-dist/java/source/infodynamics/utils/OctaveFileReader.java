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
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Array;

/**
 * Octave text file format reader.
 * 
 * Usage:
 * <ol>
 * <li>call constructor</li>
 * <li>to retrieve a variable, call the get method appropriate to the
 *    data type of that variable.</li>
 * </ol>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class OctaveFileReader {
	private String filename;
	private BufferedReader br;
	private boolean isOpen = false;
	
	static final String SCALAR_HEADER = "# type: scalar";
	static final String GLOBAL_SCALAR_HEADER = "# type: scalar";
	static final String MATRIX_HEADER = "# type: matrix";
	static final String GLOBAL_MATRIX_HEADER = "# type: global matrix";
	static final String BOOL_MATRIX_HEADER = "# type: bool matrix";
	static final String OCTAVE_DELIMITER = "[ \t]";
	
	public OctaveFileReader(String octaveFilename) throws FileNotFoundException {
		filename = octaveFilename;
		openOctaveFile();
	}
	
	private void openOctaveFile() throws FileNotFoundException {
		if (isOpen) {
			return;
		}
		br = new BufferedReader(new FileReader(filename));
		isOpen = true;
	}
	
	private void closeOctaveFile() throws IOException {
		br.close();
		isOpen = false;
	}
	
	/*
	public OctaveFileReader(FileReader fileReader) {
		br = new BufferedReader(fileReader);
	}
	
	public OctaveFileReader(File file) throws FileNotFoundException {
		br = new BufferedReader(new FileReader(file));
	}
	*/
	
	public int getInt(String name) throws Exception {
		openOctaveFile();
		try {
			for (String line = br.readLine(); line != null; line = br.readLine()) {
				if (line.equalsIgnoreCase("# name: " + name)) {
					// We have the header line for this stored variable
					// Make sure the next line indicates a scalar
					line = br.readLine();
					if (!line.equals(SCALAR_HEADER) && !line.equals(GLOBAL_SCALAR_HEADER)) {
						// The variable is not a scalar - it is not compatible with the
						//  int type
						throw new Exception("Stored Octave variable " + name + " is not of a type compatible with int");
					}
					// The next line contains the variable's value
					line = br.readLine();
					closeOctaveFile();
					// Now, since we want to allow a call to getInt even if this is a double,
					// then we'll call parse double then cast it to an int.
					double value = Double.parseDouble(line);
					return (int) value;
				}
			}
			// End of file reached and variable not found !
			throw new Exception("Stored Octave variable " + name + " not found");
		} catch (Exception e) {
			// clean up
			closeOctaveFile();
			// let the exception go
			throw e;
		}
	}
	
	public double getDouble(String name) throws Exception {
		openOctaveFile();
		try {
			for (String line = br.readLine(); line != null; line = br.readLine()) {
				if (line.equalsIgnoreCase("# name: " + name)) {
					// We have the header line for this stored variable
					// Make sure the next line indicates a scalar
					line = br.readLine();
					if (!line.equals(SCALAR_HEADER) && !line.equals(GLOBAL_SCALAR_HEADER)) {
						// The variable is not a scalar - it is not compatible with the
						//  int type
						throw new Exception("Stored Octave variable " + name + " is not of a type compatible with double");
					}
					// The next line contains the variable's value
					line = br.readLine();
					closeOctaveFile();
					return Double.parseDouble(line);
				}
			}
			// End of file reached and variable not found !
			throw new Exception("Stored Octave variable " + name + " not found");
		} catch (Exception e) {
			// clean up
			closeOctaveFile();
			throw e;
		}
	}

	public int[] getIntArray(String name) throws Exception {
		openOctaveFile();
		try {
			for (String line = br.readLine(); line != null; line = br.readLine()) {
				if (line.equalsIgnoreCase("# name: " + name)) {
					// We have the header line for this stored variable
					// Make sure the next line indicates a matrix
					line = br.readLine();
					if (!line.equals(MATRIX_HEADER) && !line.equals(GLOBAL_MATRIX_HEADER)
							&& !line.equals(BOOL_MATRIX_HEADER)) {
						// The variable is not a matrix - it is not compatible with the
						//  int [] type
						throw new Exception("Stored Octave variable " + name + " is not of a type compatible with int[] - not a matrix");
					}
					// The next line contains the number of rows
					line = br.readLine();
					String[] parts = line.split(" ");
					int rows = Integer.parseInt(parts[2]);
					// The next line contains the number of columns
					line = br.readLine();
					parts = line.split(" ");
					int columns = Integer.parseInt(parts[2]);
					// Now check that it either has one row or one column
					if (!(rows == 1) && !(columns == 1)) {
						throw new Exception("Stored Octave variable " + name + " is not of a type compatible with int[] - not a single dimension (" +
								rows + " x " + columns + ")");
					}
					int theLength = (rows > columns) ? rows : columns;
					int[] values = new int[theLength];
					// The next lines contain each row of values
					int l = 0;
					for (int r = 0; r < rows; r++) {
						// Grab the text for the next row - values separated by spaces
						line = br.readLine();
						String[] stringValues = line.split(OCTAVE_DELIMITER);
						// ignore stringValues[0] since the line starts with a delimiter
						if ((stringValues.length + 1) < columns) {
							throw new Exception("Not enough values for row " + r + " of variable " + name);
						}
						for (int c = 0; c < columns; c++) {
							// Parse the next value
							values[l] = Integer.parseInt(stringValues[c + 1]);
							l++;
						}
					}
					closeOctaveFile();
					return values;
				}
			}
			// End of file reached and variable not found !
			throw new Exception("Stored Octave variable " + name + " not found");
		} catch (Exception e) {
			// clean up
			closeOctaveFile();
			throw e;
		}
	}

	public long[] getLongArray(String name) throws Exception {
		openOctaveFile();
		try {
			for (String line = br.readLine(); line != null; line = br.readLine()) {
				if (line.equalsIgnoreCase("# name: " + name)) {
					// We have the header line for this stored variable
					// Make sure the next line indicates a matrix
					line = br.readLine();
					if (!line.equals(MATRIX_HEADER) && !line.equals(GLOBAL_MATRIX_HEADER)
							&& !line.equals(BOOL_MATRIX_HEADER)) {
						// The variable is not a matrix - it is not compatible with the
						//  int [] type
						throw new Exception("Stored Octave variable " + name + " is not of a type compatible with long[] - not a matrix");
					}
					// The next line contains the number of rows
					line = br.readLine();
					String[] parts = line.split(" ");
					int rows = Integer.parseInt(parts[2]);
					// The next line contains the number of columns
					line = br.readLine();
					parts = line.split(" ");
					int columns = Integer.parseInt(parts[2]);
					// Now check that it either has one row or one column
					if ((rows != 1) && (columns != 1)) {
						throw new Exception("Stored Octave variable " + name + " is not of a type compatible with long[] - not a single dimension (" +
								rows + " x " + columns + ")");
					}
					int theLength = (rows > columns) ? rows : columns;
					long[] values = new long[theLength];
					// The next lines contain each row of values
					int l = 0;
					for (int r = 0; r < rows; r++) {
						// Grab the text for the next row - values separated by spaces
						line = br.readLine();
						String[] stringValues = line.split(OCTAVE_DELIMITER);
						// ignore stringValues[0] since the line starts with a delimiter
						if ((stringValues.length + 1) < columns) {
							throw new Exception("Not enough values for row " + r + " of variable " + name);
						}
						for (int c = 0; c < columns; c++) {
							// Parse the next value
							values[l] = Long.parseLong(stringValues[c + 1]);
							l++;
						}
					}
					closeOctaveFile();
					return values;
				}
			}
			// End of file reached and variable not found !
			throw new Exception("Stored Octave variable " + name + " not found");
		} catch (Exception e) {
			// clean up
			closeOctaveFile();
			throw e;
		}
	}

	public double[] getDoubleArray(String name) throws Exception {
		openOctaveFile();
		try {
			for (String line = br.readLine(); line != null; line = br.readLine()) {
				if (line.equalsIgnoreCase("# name: " + name)) {
					// We have the header line for this stored variable
					// Make sure the next line indicates a matrix
					line = br.readLine();
					if (!line.equals(MATRIX_HEADER) && !line.equals(GLOBAL_MATRIX_HEADER)) {
						// The variable is not a matrix - it is not compatible with the
						//  int [] type
						throw new Exception("Stored Octave variable " + name + " is not of a type compatible with double[] - not a matrix");
					}
					// The next line contains the number of rows
					line = br.readLine();
					String[] parts = line.split(" ");
					int rows = Integer.parseInt(parts[2]);
					// The next line contains the number of columns
					line = br.readLine();
					parts = line.split(" ");
					int columns = Integer.parseInt(parts[2]);
					// Now check that it either has one row or one column
					if (!(rows == 1) && !(columns == 1)) {
						throw new Exception("Stored Octave variable " + name + " is not of a type compatible with double[] - not a single dimension (" +
								rows + " x " + columns + ")");
					}
					int theLength = (rows > columns) ? rows : columns;
					double[] values = new double[theLength];
					// The next lines contain each row of values
					int l = 0;
					for (int r = 0; r < rows; r++) {
						// Grab the text for the next row - values separated by spaces
						line = br.readLine();
						String[] stringValues = line.split(OCTAVE_DELIMITER);
						// ignore stringValues[0] since the line starts with a delimiter
						if ((stringValues.length + 1) < columns) {
							throw new Exception("Not enough values for row " + r + " of variable " + name);
						}
						for (int c = 0; c < columns; c++) {
							// Parse the next value
							values[l] = Double.parseDouble(stringValues[c + 1]);
							l++;
						}
					}
					closeOctaveFile();
					return values;
				}
			}
			// End of file reached and variable not found !
			throw new Exception("Stored Octave variable " + name + " not found");
		} catch (Exception e) {
			// clean up
			closeOctaveFile();
			throw e;
		}
	}

	public int[][] getInt2DMatrix(String name) throws Exception {
		openOctaveFile();
		try {
			for (String line = br.readLine(); line != null; line = br.readLine()) {
				if (line.equalsIgnoreCase("# name: " + name)) {
					// We have the header line for this stored variable
					// Make sure the next line indicates a matirx
					line = br.readLine();
					if (!line.equals(MATRIX_HEADER) && !line.equals(GLOBAL_MATRIX_HEADER)
							&& !line.equals(BOOL_MATRIX_HEADER)) {
						// The variable is not a matrix - it is not compatible with the
						//  int [][] type
						throw new Exception("Stored Octave variable " + name + " is not of a type compatible with int[][]");
					}
					// The next line contains the number of rows
					line = br.readLine();
					String[] parts = line.split(" ");
					int rows = Integer.parseInt(parts[2]);
					// The next line contains the number of columns
					line = br.readLine();
					parts = line.split(" ");
					int columns = Integer.parseInt(parts[2]);
					int[][] values = new int[rows][columns];
					// The next lines contain each row of values
					for (int r = 0; r < rows; r++) {
						// Grab the text for the next row - values separated by spaces
						line = br.readLine();
						String[] stringValues = line.split(OCTAVE_DELIMITER);
						int subFromLength = 0;
						if ((stringValues.length > 0) && (stringValues[0].equalsIgnoreCase(""))) {
							// Ignore stringValues[0] since the line starts with a delimiter
							subFromLength = 1;
						}
						if ((stringValues.length - subFromLength) < columns) {
							throw new Exception("Not enough values for row " + r + " of " + rows + " of variable " + name + "( actual: " +
									(stringValues.length - subFromLength) + ", required: " + columns + ")");
						}
						for (int c = 0; c < columns; c++) {
							// Parse the next value
							values[r][c] = Integer.parseInt(stringValues[c + subFromLength]);
						}
					}
					closeOctaveFile();
					return values;
				}
			}
			// End of file reached and variable not found !
			throw new Exception("Stored Octave variable " + name + " not found");
		} catch (Exception e) {
			// clean up
			closeOctaveFile();
			throw e;
		}
	}

	public double[][] getDouble2DMatrix(String name) throws Exception {
		openOctaveFile();
		try {
			for (String line = br.readLine(); line != null; line = br.readLine()) {
				if (line.equalsIgnoreCase("# name: " + name)) {
					// We have the header line for this stored variable
					// Make sure the next line indicates a matirx
					line = br.readLine();
					if (!line.equals(MATRIX_HEADER) && !line.equals(GLOBAL_MATRIX_HEADER)) {
						// The variable is not a matrix - it is not compatible with the
						//  int [][] type
						throw new Exception("Stored Octave variable " + name + " is not of a type compatible with int[][]");
					}
					// The next line contains the number of rows
					line = br.readLine();
					String[] parts = line.split(" ");
					int rows = Integer.parseInt(parts[2]);
					// The next line contains the number of columns
					line = br.readLine();
					parts = line.split(" ");
					int columns = Integer.parseInt(parts[2]);
					double[][] values = new double[rows][columns];
					// The next lines contain each row of values
					for (int r = 0; r < rows; r++) {
						// Grab the text for the next row - values separated by spaces
						line = br.readLine();
						String[] stringValues = line.split(OCTAVE_DELIMITER);
						// ignore stringValues[0] since the line starts with a delimiter
						// TODO Need to fix this in case it doesn't start with a delimiter.
						//  See the 2D integer method. Should fix other methods also. 
						if ((stringValues.length + 1) < columns) {
							throw new Exception("Not enough values for row " + r + " of variable " + name);
						}
						for (int c = 0; c < columns; c++) {
							// Parse the next value
							values[r][c] = Double.parseDouble(stringValues[c + 1]);
						}
					}
					closeOctaveFile();
					return values;
				}
			}
			// End of file reached and variable not found !
			throw new Exception("Stored Octave variable " + name + " not found");
		} catch (Exception e) {
			// clean up
			closeOctaveFile();
			throw e;
		}
	}

	public static void main(String argv[]) throws Exception {

		// Run some quick tests - ALL WORKS
		
		OctaveFileReader ofr = new OctaveFileReader("C:\\Work\\Investigations\\CATransferEntropy\\base-2\\neighbourhood-5\\rule-133976044.txt");
		// Now grab a matrix variable
		int[][] rules = ofr.getInt2DMatrix("settledExecutedRules");
		for (int r = 0; r < rules.length; r++) {
			for (int c = 0; c < rules[r].length; c++) {
				System.out.print(" " + rules[r][c]);
			}
			System.out.println();
		}
		System.out.println("rules has " + rules.length + " rows and " + rules[0].length + " columns");
		/* double[][] rulesDouble = ofr.getDouble2DMatrix("settledExecutedRules");
		for (int r = 0; r < rulesDouble.length; r++) {
			for (int c = 0; c < rulesDouble[r].length; c++) {
				System.out.print(" " + rulesDouble[r][c]);
			}
			System.out.println();
		}
		System.out.println("rulesDouble has " + rulesDouble.length + " rows and " + rulesDouble[0].length + " columns");
		*/
		System.out.println("base = " + ofr.getInt("base"));
		System.out.println("neighbourhood = " + ofr.getInt("neighbourhood"));
		
		Object obj = new int[3];
		if (obj.getClass().isArray()) {
			System.out.println("got an array of length " + Array.getLength(obj));
		}
		System.out.println(obj.getClass() + " " + obj);
	}
}
